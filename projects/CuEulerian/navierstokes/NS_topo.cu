#include "Structures.hpp"
#include "zensim/cuda/execution/ExecutionPolicy.cuh"
#include "zensim/geometry/LevelSetUtils.tpp"
#include "zensim/geometry/SparseGrid.hpp"
#include "zensim/geometry/VdbLevelSet.h"
#include "zensim/omp/execution/ExecutionPolicy.hpp"
#include "zensim/profile/CppTimers.hpp"
#include "zensim/zpc_tpls/fmt/color.h"
#include "zensim/zpc_tpls/fmt/format.h"

#include <zeno/types/ListObject.h>
#include <zeno/types/NumericObject.h>
#include <zeno/types/PrimitiveObject.h>

#include <zeno/VDBGrid.h>

#include "../scheme.hpp"
#include "../utils.cuh"

namespace zeno {

struct ZSExtendSparseGrid : INode {

    template <typename PredT> void refit(ZenoSparseGrid *nsgridPtr, zs::SmallString tag, PredT pred) {
        using namespace zs;
        static constexpr auto space = execspace_e::cuda;
        namespace cg = ::cooperative_groups;
        auto pol = cuda_exec();
        auto &spg = nsgridPtr->getSparseGrid();
        // make sure spg.block_size % 32 == 0

        auto nbs = spg.numBlocks();
        using Ti = RM_CVREF_T(nbs);

        Vector<Ti> marks{spg.get_allocator(), nbs + 1}, offsets{spg.get_allocator(), nbs + 1};
        marks.reset(0);

        pol(range(nbs * 32), [spg = proxy<space>(spg), tagOffset = spg.getPropertyOffset(tag),
                              marks = proxy<space>(marks), pred] __device__(std::size_t i) mutable {
            auto tile = cg::tiled_partition<32>(cg::this_thread_block());
            auto bno = i / 32;
            auto cellno = tile.thread_rank();

            while (cellno < spg.block_size) {
                if (tile.ballot(pred(spg(tagOffset, bno, cellno))))
                    break;
                cellno += 32;
            }
            if (tile.thread_rank() == 0 && cellno < spg.block_size)
                marks[bno] = 1;
        });

        exclusive_scan(pol, std::begin(marks), std::end(marks), std::begin(offsets));
        auto newNbs = offsets.getVal(nbs);
        fmt::print("compacting {} blocks to {} active blocks.\n", nbs, newNbs);

        /// @brief compact active blocks
        // grid
        auto &grid = spg._grid;
        auto dstgrid = grid;
        // table
        auto &table = spg._table;
        table.reset(false);
        table._cnt.setVal(newNbs);
        auto &keys = table._activeKeys;
        auto newKeys = keys;

        pol(range(nbs * spg._table.bucket_size),
            [grid = proxy<space>(grid), dstgrid = proxy<space>(dstgrid), marks = proxy<space>(marks),
             newKeys = proxy<space>(newKeys), keys = proxy<space>(keys), table = proxy<space>(table),
             offsets = proxy<space>(offsets),
             bs_c = wrapv<RM_CVREF_T(spg)::block_size>{}] __device__(std::size_t i) mutable {
                constexpr auto block_size = decltype(bs_c)::value;
                constexpr auto bucket_size = RM_CVREF_T(table)::bucket_size;
                static_assert(block_size % bucket_size == 0, "block_size should be a multiple of bucket_size");
                auto tile = cg::tiled_partition<bucket_size>(cg::this_thread_block());
                auto bno = i / bucket_size;
                if (marks[bno] == 0)
                    return;
                auto dstBno = offsets[bno];
                // grid
                for (auto cellno = tile.thread_rank(); cellno < block_size; cellno += bucket_size) {
                    for (int chn = 0; chn != grid.numChannels(); ++chn)
                        dstgrid(chn, dstBno, cellno) = grid(chn, bno, cellno);
                }
                // table
                auto bcoord = keys[bno];
                table.tile_insert(tile, bcoord, dstBno, false); // do not enqueue key, hence set false
                if (tile.thread_rank() == 0)
                    newKeys[dstBno] = bcoord;
            });
        grid = std::move(dstgrid);
        keys = std::move(newKeys);
    }

    template <typename PredT>
    void extend(ZenoSparseGrid *nsgridPtr, zs::SmallString tag, std::size_t &nbsOffset, PredT pred) {
        using namespace zs;
        static constexpr auto space = execspace_e::cuda;
        namespace cg = ::cooperative_groups;
        auto pol = cuda_exec();
        auto &spg = nsgridPtr->getSparseGrid();
        // make sure spg.block_size % 32 == 0

        auto nbs = spg.numBlocks() - nbsOffset;
        if (nbs == 0)
            return;
        // worst case is that all candidate blocks activate all surrounding neighbors
        spg.resize(pol, nbs * 26 + nbsOffset);

        // zeno::log_info("currently {} blocks (offset {}), resizing to {}\n", nbsOffset + nbs, nbsOffset,
        //                nbs * 26 + nbsOffset);

        if (!spg._grid.hasProperty(tag))
            throw std::runtime_error(fmt::format("property [{}] not exist!", tag.asString()));

        pol(range(nbs * spg._table.bucket_size), [spg = proxy<space>(spg), tagOffset = spg.getPropertyOffset(tag),
                                                  nbsOffset, pred] __device__(std::size_t i) mutable {
            auto tile = cg::tiled_partition<RM_CVREF_T(spg._table)::bucket_size>(cg::this_thread_block());
            auto bno = i / spg._table.bucket_size + nbsOffset;
            auto cellno = tile.thread_rank();
            // searching for active voxels within this block

            while (cellno < spg.block_size) {
                if (tile.ballot(pred(spg(tagOffset, bno, cellno))))
                    break;
                cellno += spg._table.bucket_size;
            }
            if (cellno < spg.block_size) {
                auto bcoord = spg.iCoord(bno, 0);
                for (auto loc : ndrange<3>(3)) {
                    auto dir = make_vec<int>(loc) - 1;
                    // spg._table.insert(bcoord + dir * spg.side_length);
                    spg._table.tile_insert(tile, bcoord + dir * spg.side_length, RM_CVREF_T(spg._table)::sentinel_v,
                                           true);
                }
            }
        });
        // [ nbsOffset | nbsOffset + nbs | spg.numBlocks() ]
        nbsOffset += nbs;
        auto newNbs = spg.numBlocks();
        newNbs -= nbsOffset;
        if (newNbs > 0)
            zs::memset(mem_device, (void *)spg._grid.tileOffset(nbsOffset), 0,
                       (std::size_t)newNbs * spg._grid.tileBytes());

        if (tag == "sdf")
            pol(range(newNbs * spg.block_size),
                [dx = spg.voxelSize()[0], spg = proxy<space>(spg), sdfOffset = spg.getPropertyOffset("sdf"),
                 blockOffset = nbsOffset * spg.block_size] __device__(std::size_t cellno) mutable {
                    spg(sdfOffset, blockOffset + cellno) = 3 * dx;
                });
    }

    void apply() override {
        auto zsSPG = get_input<ZenoSparseGrid>("NSGrid");
        auto tag = get_input2<std::string>("Attribute");
        auto nlayers = get_input2<int>("layers");
        auto needRefit = get_input2<bool>("refit");

        std::size_t nbs = 0;
        int opt = 0;
        if (tag == "rho")
            opt = 1;
        else if (tag == "sdf")
            opt = 2;

        if (needRefit && opt != 0) {
            if (opt == 1)
                refit(zsSPG.get(), src_tag(zsSPG, tag),
                      [] __device__(float v) -> bool { return v > zs::limits<float>::epsilon() * 10; });
            else if (opt == 2)
                refit(zsSPG.get(), src_tag(zsSPG, tag),
                      [dx = zsSPG->getSparseGrid().voxelSize()[0]] __device__(float v) -> bool { return v < 2 * dx; });
            opt = 0;
        }

        while (nlayers-- > 0) {
            if (opt == 0)
                extend(zsSPG.get(), src_tag(zsSPG, tag), nbs, [] __device__(float v) { return true; });
            else if (opt == 1)
                extend(zsSPG.get(), src_tag(zsSPG, tag), nbs,
                       [] __device__(float v) -> bool { return v > zs::limits<float>::epsilon() * 10; });
            else if (opt == 2)
                extend(zsSPG.get(), src_tag(zsSPG, tag), nbs,
                       [dx = zsSPG->getSparseGrid().voxelSize()[0]] __device__(float v) -> bool { return v < 2 * dx; });
            opt = 0; // always active since
        }

        set_output("NSGrid", zsSPG);
    }
};

ZENDEFNODE(ZSExtendSparseGrid,
           {/* inputs: */
            {"NSGrid", {"enum rho sdf", "Attribute", "rho"}, {"bool", "refit", "1"}, {"int", "layers", "2"}},
            /* outputs: */
            {"NSGrid"},
            /* params: */
            {},
            /* category: */
            {"Eulerian"}});

struct ZSMaintainSparseGrid : INode {
    template <typename PredT> void maintain(ZenoSparseGrid *nsgridPtr, zs::SmallString tag, PredT pred, int nlayers) {
        using namespace zs;
        static constexpr auto space = execspace_e::cuda;
        namespace cg = ::cooperative_groups;
        auto pol = cuda_exec();
        auto &spg = nsgridPtr->getSparseGrid();

        if (!spg._grid.hasProperty(tag))
            throw std::runtime_error(fmt::format("property [{}] not exist!", tag.asString()));

        auto nbs = spg.numBlocks();
        using Ti = RM_CVREF_T(nbs);

        Vector<Ti> marks{spg.get_allocator(), nbs + 1}, offsets{spg.get_allocator(), nbs + 1};
        marks.reset(0);

        static_assert(RM_CVREF_T(spg)::block_size % 32 == 0, "block size should be a multiple of 32.");

        /// @brief mark active block entries
        pol(range(nbs * 32), [spg = proxy<space>(spg), tagOffset = spg.getPropertyOffset(tag),
                              marks = proxy<space>(marks), pred] __device__(std::size_t i) mutable {
            auto tile = cg::tiled_partition<32>(cg::this_thread_block());
            auto bno = i / 32;
            auto cellno = tile.thread_rank();

            while (cellno < spg.block_size) {
                if (tile.ballot(pred(spg(tagOffset, bno, cellno))))
                    break;
                cellno += 32;
            }
            if (tile.thread_rank() == 0 && cellno < spg.block_size)
                marks[bno] = 1;
        });

        exclusive_scan(pol, std::begin(marks), std::end(marks), std::begin(offsets));
        auto newNbs = offsets.getVal(nbs);
        fmt::print("compacting {} blocks to {} active blocks.\n", nbs, newNbs);

        /// @brief compact active block entries
        // table
        auto &table = spg._table;
        table.reset(false);
        table._cnt.setVal(newNbs);
        // backup previous block entries, nbs is the previous count of blocks
        auto prevKeys = table._activeKeys;
        auto &keys = table._activeKeys;

        pol(range(nbs * spg._table.bucket_size),
            [marks = proxy<space>(marks), newKeys = proxy<space>(keys), keys = proxy<space>(prevKeys),
             table = proxy<space>(table), offsets = proxy<space>(offsets),
             bs_c = wrapv<RM_CVREF_T(spg)::block_size>{}] __device__(std::size_t i) mutable {
                constexpr auto block_size = decltype(bs_c)::value;
                constexpr auto bucket_size = RM_CVREF_T(table)::bucket_size;
                static_assert(block_size % bucket_size == 0, "block_size should be a multiple of bucket_size");
                auto tile = cg::tiled_partition<bucket_size>(cg::this_thread_block());
                auto bno = i / bucket_size;
                if (marks[bno] == 0)
                    return;
                auto dstBno = offsets[bno];
                // table
                auto bcoord = keys[bno];
                table.tile_insert(tile, bcoord, dstBno, false); // do not enqueue key, hence set false
                if (tile.thread_rank() == 0)
                    newKeys[dstBno] = bcoord;
            });

        // grid
        /// @note backup the grid ahead
        auto &grid = spg._grid;
        auto prevGrid = grid;

        /// @brief iteratively expand the active domain
        Ti nbsOffset = 0;
        while (nlayers-- > 0 && newNbs > 0) {
            // reserve enough memory for expanded grid and table
            spg.resize(pol, newNbs * 7 + nbsOffset);
            // extend one layer
            pol(range(newNbs * spg._table.bucket_size),
                [spg = proxy<space>(spg), tagOffset = spg.getPropertyOffset(tag),
                 nbsOffset] __device__(std::size_t i) mutable {
                    auto tile = cg::tiled_partition<RM_CVREF_T(spg._table)::bucket_size>(cg::this_thread_block());
                    auto bno = i / spg._table.bucket_size + nbsOffset;
                    auto bcoord = spg.iCoord(bno, 0);
                    {
                        for (int d = 0; d != 3; ++d) {
                            auto dir = zs::vec<int, 3>::zeros();
                            dir[d] = -1;
                            spg._table.tile_insert(tile, bcoord + dir * spg.side_length, RM_CVREF_T(spg._table)::sentinel_v,
                                               true);
                            dir[d] = 1;
                            spg._table.tile_insert(tile, bcoord + dir * spg.side_length, RM_CVREF_T(spg._table)::sentinel_v,
                                               true);
                        }
                    }
                    #if 0
                    for (auto loc : ndrange<3>(3)) {
                        auto dir = make_vec<int>(loc) - 1;
                        // spg._table.insert(bcoord + dir * spg.side_length);
                        spg._table.tile_insert(tile, bcoord + dir * spg.side_length, RM_CVREF_T(spg._table)::sentinel_v,
                                               true);
                    }
                    #endif
                });
            // slide the window
            nbsOffset += newNbs;
            newNbs = spg.numBlocks() - nbsOffset;

            // initialize newly added blocks
            if (newNbs > 0) {
                zs::memset(mem_device, (void *)spg._grid.tileOffset(nbsOffset), 0,
                           (std::size_t)newNbs * spg._grid.tileBytes());

                if (tag == "sdf") {
                    // special treatment for "sdf" property
                    pol(range(newNbs * spg.block_size),
                        [dx = spg.voxelSize()[0], spg = proxy<space>(spg), sdfOffset = spg.getPropertyOffset("sdf"),
                         blockOffset = nbsOffset * spg.block_size] __device__(std::size_t cellno) mutable {
                            spg(sdfOffset, blockOffset + cellno) = 3 * dx;
                        });
                }
            }
        }

        /// @brief relocate original grid data to the new sparse grid
        pol(range(nbs * spg._table.bucket_size), [grid = proxy<space>(prevGrid), spg = proxy<space>(spg),
                                                  keys = proxy<space>(prevKeys)] __device__(std::size_t i) mutable {
            constexpr auto bucket_size = RM_CVREF_T(table)::bucket_size;
            auto tile = cg::tiled_partition<bucket_size>(cg::this_thread_block());
            auto bno = i / bucket_size;
            auto bcoord = keys[bno];
            auto dstBno = spg._table.tile_query(tile, bcoord);
            if (dstBno == spg._table.sentinel_v)
                return;
            // table
            for (auto cellno = tile.thread_rank(); cellno < spg.block_size; cellno += bucket_size) {
                for (int chn = 0; chn != grid.numChannels(); ++chn)
                    spg._grid(chn, dstBno, cellno) = grid(chn, bno, cellno);
            }
        });

        /// @brief adjust multigrid accordingly
        // grid
        nbs = spg.numBlocks();
        auto &spg1 = nsgridPtr->spg1;
        spg1.resize(pol, nbs);
        auto &spg2 = nsgridPtr->spg2;
        spg2.resize(pol, nbs);
        auto &spg3 = nsgridPtr->spg3;
        spg3.resize(pol, nbs);
        // table
        {
            const auto &table = spg._table;
            auto &table1 = spg1._table;
            auto &table2 = spg2._table;
            auto &table3 = spg3._table;
            table1.reset(true);
            table1._cnt.setVal(nbs);
            table2.reset(true);
            table2._cnt.setVal(nbs);
            table3.reset(true);
            table3._cnt.setVal(nbs);
            pol(range(nbs), [table = proxy<space>(table), tab1 = proxy<space>(table1), tab2 = proxy<space>(table2),
                             tab3 = proxy<space>(table3)] __device__(std::size_t i) mutable {
                auto bcoord = table._activeKeys[i];
                tab1.insertUnsafe(bcoord / 2, i, true);
                tab2.insertUnsafe(bcoord / 4, i, true);
                tab3.insertUnsafe(bcoord / 8, i, true);
            });
        }
    }

    void apply() override {
        auto zsSPG = get_input<ZenoSparseGrid>("NSGrid");
        auto tag = get_input2<std::string>("Attribute");
        auto nlayers = get_input2<int>("layers");
        auto needRefit = get_input2<bool>("refit");

        int opt = 0;
        if (needRefit) {
            if (tag == "rho")
                opt = 1;
            else if (tag == "sdf")
                opt = 2;
        }

        if (opt == 0)
            maintain(
                zsSPG.get(), src_tag(zsSPG, tag), [] __device__(float v) { return true; }, nlayers);
        else if (opt == 1)
            maintain(
                zsSPG.get(), src_tag(zsSPG, tag),
                [] __device__(float v) -> bool { return v > zs::limits<float>::epsilon() * 10; }, nlayers);
        else if (opt == 2)
            maintain(
                zsSPG.get(), src_tag(zsSPG, tag),
                [dx = zsSPG->getSparseGrid().voxelSize()[0]] __device__(float v) -> bool { return v < 2 * dx; },
                nlayers);

        set_output("NSGrid", zsSPG);
    }
};

ZENDEFNODE(ZSMaintainSparseGrid,
           {/* inputs: */
            {"NSGrid", {"enum rho sdf", "Attribute", "rho"}, {"bool", "refit", "1"}, {"int", "layers", "2"}},
            /* outputs: */
            {"NSGrid"},
            /* params: */
            {},
            /* category: */
            {"Eulerian"}});

} // namespace zeno