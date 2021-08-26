#include <openvdb/openvdb.h>
#include <openvdb/tools/Dense.h>
#include <fdb/openvdb.h>
#include <fdb/vec.h>

namespace fdb {

template <class GridT, class ValT>
void impl_write_dense_vdb
        ( std::string_view path
        , std::function<ValT(vec3I)> sampler
        , vec3I size
        )
{
    openvdb::tools::Dense<typename GridT::ValueType> dens(openvdb::Coord(size[0], size[1], size[2]));
    for (uint32_t z = 0; z < size[2]; z++) {
        for (uint32_t y = 0; y < size[1]; y++) {
            for (uint32_t x = 0; x < size[0]; x++) {
                auto val = sampler({x, y, z});
                dens.setValue(x, y, z, vec_to_other<typename GridT::ValueType>(val));
            }
        }
    }
    auto grid = GridT::create();
    typename GridT::ValueType tolerance{0};
    openvdb::tools::copyFromDense(dens, grid->tree(), tolerance);
    openvdb::io::File((std::string)path).write({grid});
}

void write_dense_vdb
    ( std::string_view path
    , std::function<float(vec3I)> sampler
    , vec3I size
    )
{
    return impl_write_dense_vdb<openvdb::FloatGrid>(path, sampler, size);
}

void write_dense_vdb
    ( std::string_view path
    , std::function<vec3f(vec3I)> sampler
    , vec3I size
    )
{
    return impl_write_dense_vdb<openvdb::Vec3fGrid>(path, sampler, size);
}

}
