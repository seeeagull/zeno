// https://github.com/alembic/alembic/blob/master/lib/Alembic/AbcGeom/Tests/PolyMeshTest.cpp
// WHY THE FKING ALEMBIC OFFICIAL GIVES NO DOC BUT ONLY "TESTS" FOR ME TO LEARN THEIR FKING LIB
#include <zeno/zeno.h>
#include <zeno/utils/logger.h>
#include <zeno/types/StringObject.h>
#include <zeno/types/PrimitiveObject.h>
#include <zeno/types/PrimitiveTools.h>
#include <zeno/types/NumericObject.h>
#include <zeno/types/UserData.h>
#include <zeno/extra/GlobalState.h>
#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreAbstract/All.h>
#include <Alembic/AbcCoreOgawa/All.h>
#include <Alembic/AbcCoreHDF5/All.h>
#include <Alembic/Abc/ErrorHandler.h>
#include "ABCTree.h"
#include "zeno/types/DictObject.h"
#include <cstring>
#include <cstdio>
#include <filesystem>
#include <zeno/utils/string.h>
#include <zeno/utils/scope_exit.h>

#ifdef ZENO_WITH_PYTHON3
    #include <Python.h>
#endif

using namespace Alembic::AbcGeom;

namespace zeno {
static int clamp(int i, int _min, int _max) {
    if (i < _min) {
        return _min;
    } else if (i > _max) {
        return _max;
    } else {
        return i;
    }
}

static void set_time_info(UserData &ud, TimeSamplingType tst, float start, int sample_count) {
    float time_per_cycle = tst.getTimePerCycle();
    if (tst.isUniform()) {
        ud.set2("_abc_time_sampling_type", "Uniform");
    }
    else if (tst.isCyclic()) {
        ud.set2("_abc_time_sampling_type", "Cyclic");
    }
    else if (tst.isAcyclic()) {
        ud.set2("_abc_time_sampling_type", "Acyclic");
    }
    ud.set2("_abc_start_time", float(start));
    ud.set2("_abc_sample_count", sample_count);
    ud.set2("_abc_time_per_cycle", time_per_cycle);
    if (time_per_cycle > 0) {
        ud.set2("_abc_time_fps", 1.0f / time_per_cycle);
    }
    else {
        ud.set2("_abc_time_fps", 0.0f);
    }
}
static void read_velocity(std::shared_ptr<PrimitiveObject> prim, V3fArraySamplePtr marr, bool read_done) {
    if (marr == nullptr) {
        return;
    }
    if (marr->size() > 0) {
        if (!read_done) {
            log_info("[alembic] totally {} velocities", marr->size());
        }
        auto &parr = prim->add_attr<vec3f>("v");
        for (size_t i = 0; i < marr->size(); i++) {
            auto const &val = (*marr)[i];
            parr[i] = {val[0], val[1], val[2]};
        }
    }
}

static void read_attributes(std::shared_ptr<PrimitiveObject> prim, ICompoundProperty arbattrs, const ISampleSelector &iSS, bool read_done) {
    if (!arbattrs) {
        return;
    }
    size_t numProps = arbattrs.getNumProperties();
    for (auto i = 0; i < numProps; i++) {
        PropertyHeader p = arbattrs.getPropertyHeader(i);
        if (IFloatGeomParam::matches(p)) {
            IFloatGeomParam param(arbattrs, p.getName());

            IFloatGeomParam::Sample samp = param.getIndexedValue(iSS);
            std::vector<float> data;
            data.resize(samp.getVals()->size());
            for (auto i = 0; i < samp.getVals()->size(); i++) {
                data[i] = samp.getVals()->get()[i];
            }
            if (!read_done) {
                log_info("[alembic] float attr {}, len {}.", p.getName(), data.size());
            }

            if (prim->verts.size() == data.size()) {
                auto &attr = prim->add_attr<float>(p.getName());
                for (auto i = 0; i < prim->verts.size(); i++) {
                    attr[i] = data[i];
                }
            }
            else if (prim->verts.size() * 3 == data.size()) {
                auto &attr = prim->add_attr<zeno::vec3f>(p.getName());
                for (auto i = 0; i < prim->verts.size(); i++) {
                    attr[i] = { data[ 3 * i], data[3 * i + 1], data[3 * i + 2]};
                }
            }
            else if (prim->polys.size() == data.size()) {
                auto &attr = prim->polys.add_attr<float>(p.getName());
                for (auto i = 0; i < prim->polys.size(); i++) {
                    attr[i] = data[i];
                }
            }
            else if (prim->polys.size() * 3 == data.size()) {
                auto &attr = prim->polys.add_attr<zeno::vec3f>(p.getName());
                for (auto i = 0; i < prim->polys.size(); i++) {
                    attr[i] = { data[ 3 * i], data[3 * i + 1], data[3 * i + 2]};
                }
            }
            else if (prim->loops.size() == data.size()) {
                auto &attr = prim->loops.add_attr<float>(p.getName());
                for (auto i = 0; i < prim->loops.size(); i++) {
                    attr[i] = data[i];
                }
            }
            else if (prim->loops.size() * 3 == data.size()) {
                auto &attr = prim->loops.add_attr<zeno::vec3f>(p.getName());
                for (auto i = 0; i < prim->loops.size(); i++) {
                    attr[i] = { data[ 3 * i], data[3 * i + 1], data[3 * i + 2]};
                }
            }
            else {
                if (!read_done) {
                    log_warn("[alembic] can not load float attr {}: {}.", p.getName(), data.size());
                }
            }
        }
        else if (IInt32GeomParam::matches(p)) {
            IInt32GeomParam param(arbattrs, p.getName());

            IInt32GeomParam::Sample samp = param.getIndexedValue(iSS);
            std::vector<int> data;
            data.resize(samp.getVals()->size());
            for (auto i = 0; i < samp.getVals()->size(); i++) {
                data[i] = samp.getVals()->get()[i];
            }
            if (!read_done) {
                log_info("[alembic] i32 attr {}, len {}.", p.getName(), data.size());
            }

            if (prim->verts.size() == data.size()) {
                auto &attr = prim->add_attr<int>(p.getName());
                for (auto i = 0; i < prim->verts.size(); i++) {
                    attr[i] = data[i];
                }
            }
            else if (prim->loops.size() == data.size()) {
                auto &attr = prim->loops.add_attr<int>(p.getName());
                for (auto i = 0; i < prim->loops.size(); i++) {
                    attr[i] = data[i];
                }
            }
            else if (prim->polys.size() == data.size()) {
                auto &attr = prim->polys.add_attr<int>(p.getName());
                for (auto i = 0; i < prim->polys.size(); i++) {
                    attr[i] = data[i];
                }
            }
            else {
                if (!read_done) {
                    log_warn("[alembic] can not load int attr {}:{}.", p.getName(), data.size());
                }
            }
        }
        else if (IV3fGeomParam::matches(p)) {
            IV3fGeomParam param(arbattrs, p.getName());
            if (!read_done) {
                log_info("[alembic] vec3f attr {}.", p.getName());
            }
            IV3fGeomParam::Sample samp = param.getIndexedValue(iSS);
            if (prim->verts.size() == samp.getVals()->size()) {
                auto &attr = prim->add_attr<zeno::vec3f>(p.getName());
                for (auto i = 0; i < prim->verts.size(); i++) {
                    auto v = samp.getVals()->get()[i];
                    attr[i] = {v[0], v[1], v[2]};
                }
            }
            else if (prim->loops.size() == samp.getVals()->size()) {
                auto &attr = prim->loops.add_attr<zeno::vec3f>(p.getName());
                for (auto i = 0; i < prim->loops.size(); i++) {
                    auto v = samp.getVals()->get()[i];
                    attr[i] = {v[0], v[1], v[2]};
                }
            }
            else if (prim->polys.size() == samp.getVals()->size()) {
                auto &attr = prim->polys.add_attr<zeno::vec3f>(p.getName());
                for (auto i = 0; i < prim->polys.size(); i++) {
                    auto v = samp.getVals()->get()[i];
                    attr[i] = {v[0], v[1], v[2]};
                }
            }
            else {
                if (!read_done) {
                    log_warn("[alembic] can not load vec3f attr {}:{}.", p.getName(), int(samp.getVals()->size()));
                }
            }
        }
        else if (IN3fGeomParam::matches(p)) {
            if (!read_done) {
                log_info("[alembic] IN3fGeomParam attr {}.", p.getName());
            }
            IN3fGeomParam param(arbattrs, p.getName());
            IN3fGeomParam::Sample samp = param.getIndexedValue(iSS);
            if (prim->verts.size() == samp.getVals()->size()) {
                auto &attr = prim->add_attr<zeno::vec3f>(p.getName());
                for (auto i = 0; i < prim->verts.size(); i++) {
                    auto v = samp.getVals()->get()[i];
                    attr[i] = {v[0], v[1], v[2]};
                }
            }
            else if (prim->loops.size() == samp.getVals()->size()) {
                auto &attr = prim->loops.add_attr<zeno::vec3f>(p.getName());
                for (auto i = 0; i < prim->loops.size(); i++) {
                    auto v = samp.getVals()->get()[i];
                    attr[i] = {v[0], v[1], v[2]};
                }
            }
            else if (prim->polys.size() == samp.getVals()->size()) {
                auto &attr = prim->polys.add_attr<zeno::vec3f>(p.getName());
                for (auto i = 0; i < prim->polys.size(); i++) {
                    auto v = samp.getVals()->get()[i];
                    attr[i] = {v[0], v[1], v[2]};
                }
            }
            else {
                if (!read_done) {
                    log_warn("[alembic] can not load N3f attr {}:{}.", p.getName(), int(samp.getVals()->size()));
                }
            }
        }
        else if (IC3fGeomParam::matches(p)) {
            if (!read_done) {
                log_info("[alembic] IC3fGeomParam attr {}.", p.getName());
            }
            IC3fGeomParam param(arbattrs, p.getName());
            IC3fGeomParam::Sample samp = param.getIndexedValue(iSS);
            if (prim->verts.size() == samp.getVals()->size()) {
                auto &attr = prim->add_attr<zeno::vec3f>(p.getName());
                for (auto i = 0; i < prim->verts.size(); i++) {
                    auto v = samp.getVals()->get()[i];
                    attr[i] = {v[0], v[1], v[2]};
                }
            }
            else if (prim->loops.size() == samp.getVals()->size()) {
                auto &attr = prim->loops.add_attr<zeno::vec3f>(p.getName());
                for (auto i = 0; i < prim->loops.size(); i++) {
                    auto v = samp.getVals()->get()[i];
                    attr[i] = {v[0], v[1], v[2]};
                }
            }
            else if (prim->polys.size() == samp.getVals()->size()) {
                auto &attr = prim->polys.add_attr<zeno::vec3f>(p.getName());
                for (auto i = 0; i < prim->polys.size(); i++) {
                    auto v = samp.getVals()->get()[i];
                    attr[i] = {v[0], v[1], v[2]};
                }
            }
            else {
                if (!read_done) {
                    log_warn("[alembic] can not load C3f attr {}:{}.", p.getName(), int(samp.getVals()->size()));
                }
            }
        }
        else {
            if (!read_done) {
                log_warn("[alembic] can not load attr {}..", p.getName());
            }
        }
    }
}

static void read_user_data(std::shared_ptr<PrimitiveObject> prim, ICompoundProperty arbattrs, const ISampleSelector &iSS, bool read_done) {
    if (!arbattrs) {
        return;
    }
    size_t numProps = arbattrs.getNumProperties();
    for (auto i = 0; i < numProps; i++) {
        PropertyHeader p = arbattrs.getPropertyHeader(i);
        if (IFloatProperty::matches(p)) {
            IFloatProperty param(arbattrs, p.getName());

            float v = param.getValue(iSS);
            prim->userData().set2(p.getName(), v);
        }
        else if (IInt32Property::matches(p)) {
            IInt32Property param(arbattrs, p.getName());

            int v = param.getValue(iSS);
            prim->userData().set2(p.getName(), v);
        }
        else if (IV2fProperty::matches(p)) {
            IV2fProperty param(arbattrs, p.getName());

            auto v = param.getValue(iSS);
            prim->userData().set2(p.getName(), vec2f(v[0], v[1]));
        }
        else if (IV3fProperty::matches(p)) {
            IV3fProperty param(arbattrs, p.getName());

            auto v = param.getValue(iSS);
            prim->userData().set2(p.getName(), vec3f(v[0], v[1], v[2]));
        }
        else if (IV2iProperty::matches(p)) {
            IV2iProperty param(arbattrs, p.getName());

            auto v = param.getValue(iSS);
            prim->userData().set2(p.getName(), vec2i(v[0], v[1]));
        }
        else if (IV3iProperty::matches(p)) {
            IV3iProperty param(arbattrs, p.getName());

            auto v = param.getValue(iSS);
            prim->userData().set2(p.getName(), vec3i(v[0], v[1], v[2]));
        }
        else if (IStringProperty::matches(p)) {
            IStringProperty param(arbattrs, p.getName());

            auto value = param.getValue(iSS);
            prim->userData().set2(p.getName(), value);
        }
        else if (IBoolProperty::matches(p)) {
            IBoolProperty param(arbattrs, p.getName());

            auto value = param.getValue(iSS);
            prim->userData().set2(p.getName(), int(value));
        }
        else if (IInt16Property::matches(p)) {
            IInt16Property param(arbattrs, p.getName());

            auto value = param.getValue(iSS);
            prim->userData().set2(p.getName(), int(value));
        }
        else {
            if (!read_done) {
                log_warn("[alembic] can not load user data {}..", p.getName());
            }
        }
    }
}

static std::shared_ptr<PrimitiveObject> foundABCMesh(Alembic::AbcGeom::IPolyMeshSchema &mesh, int frameid, bool read_done, bool read_face_set) {
    auto prim = std::make_shared<PrimitiveObject>();

    std::shared_ptr<Alembic::AbcCoreAbstract::v12::TimeSampling> time = mesh.getTimeSampling();
    float time_per_cycle =  time->getTimeSamplingType().getTimePerCycle();
    double start = time->getStoredTimes().front();
    int start_frame = std::lround(start / time_per_cycle );
    set_time_info(prim->userData(), time->getTimeSamplingType(), start, int(mesh.getNumSamples()));

    int sample_index = clamp(frameid - start_frame, 0, (int)mesh.getNumSamples() - 1);
    ISampleSelector iSS = Alembic::Abc::v12::ISampleSelector((Alembic::AbcCoreAbstract::index_t)sample_index);
    Alembic::AbcGeom::IPolyMeshSchema::Sample mesamp = mesh.getValue(iSS);

    if (auto marr = mesamp.getPositions()) {
        if (!read_done) {
            log_debug("[alembic] totally {} positions", marr->size());
        }
        auto &parr = prim->verts;
        for (size_t i = 0; i < marr->size(); i++) {
            auto const &val = (*marr)[i];
            parr.emplace_back(val[0], val[1], val[2]);
        }
    }

    read_velocity(prim, mesamp.getVelocities(), read_done);
    if (auto nrm = mesh.getNormalsParam()) {
        auto nrmsamp =
                nrm.getIndexedValue(Alembic::Abc::v12::ISampleSelector((Alembic::AbcCoreAbstract::index_t)sample_index));
        int value_size = (int)nrmsamp.getVals()->size();
        if (value_size == prim->verts.size()) {
            auto &nrms = prim->verts.add_attr<vec3f>("nrm");
            auto marr = nrmsamp.getVals();
            for (size_t i = 0; i < marr->size(); i++) {
                auto const &n = (*marr)[i];
                nrms[i] = {n[0], n[1], n[2]};
            }
        }
    }

    if (auto marr = mesamp.getFaceIndices()) {
        if (!read_done) {
            log_debug("[alembic] totally {} face indices", marr->size());
        }
        auto &parr = prim->loops;
        for (size_t i = 0; i < marr->size(); i++) {
            int ind = (*marr)[i];
            parr.push_back(ind);
        }
    }

    if (auto marr = mesamp.getFaceCounts()) {
        if (!read_done) {
            log_debug("[alembic] totally {} faces", marr->size());
        }
        auto &loops = prim->loops;
        auto &parr = prim->polys;
        int base = 0;
        for (size_t i = 0; i < marr->size(); i++) {
            int cnt = (*marr)[i];
            parr.emplace_back(base, cnt);
            base += cnt;
        }
    }
    if (auto uv = mesh.getUVsParam()) {
        auto uvsamp =
            uv.getIndexedValue(Alembic::Abc::v12::ISampleSelector((Alembic::AbcCoreAbstract::index_t)sample_index));
        int value_size = (int)uvsamp.getVals()->size();
        int index_size = (int)uvsamp.getIndices()->size();
        if (!read_done) {
            log_debug("[alembic] totally {} uv value", value_size);
            log_debug("[alembic] totally {} uv indices", index_size);
            if (prim->loops.size() == index_size) {
                log_debug("[alembic] uv per face");
            } else if (prim->verts.size() == index_size) {
                log_debug("[alembic] uv per vertex");
            } else {
                log_error("[alembic] error uv indices");
            }
        }
        prim->uvs.resize(value_size);
        {
            auto marr = uvsamp.getVals();
            for (size_t i = 0; i < marr->size(); i++) {
                auto const &val = (*marr)[i];
                prim->uvs[i] = {val[0], val[1]};
            }
        }
        if (prim->loops.size() == index_size) {
            prim->loops.add_attr<int>("uvs");
            for (auto i = 0; i < prim->loops.size(); i++) {
                prim->loops.attr<int>("uvs")[i] = (*uvsamp.getIndices())[i];
            }
        }
        else if (prim->verts.size() == index_size) {
            prim->loops.add_attr<int>("uvs");
            for (auto i = 0; i < prim->loops.size(); i++) {
                prim->loops.attr<int>("uvs")[i] = prim->loops[i];
            }
        }
    }
    if (!prim->loops.has_attr("uvs")) {
        if (!read_done) {
            log_warn("[alembic] Not found uv, auto fill zero.");
        }
        prim->uvs.resize(1);
        prim->uvs[0] = zeno::vec2f(0, 0);
        prim->loops.add_attr<int>("uvs");
        for (auto i = 0; i < prim->loops.size(); i++) {
            prim->loops.attr<int>("uvs")[i] = 0;
        }
    }
    ICompoundProperty arbattrs = mesh.getArbGeomParams();
    read_attributes(prim, arbattrs, iSS, read_done);
    ICompoundProperty usrData = mesh.getUserProperties();
    read_user_data(prim, usrData, iSS, read_done);

    if (read_face_set) {
        auto &faceset = prim->polys.add_attr<int>("faceset");
        std::fill(faceset.begin(), faceset.end(), -1);
        auto &ud = prim->userData();
        std::vector<std::string> faceSetNames;
        mesh.getFaceSetNames(faceSetNames);
        ud.set2("faceset_count", int(faceSetNames.size()));
        for (auto i = 0; i < faceSetNames.size(); i++) {
            auto n = faceSetNames[i];
            ud.set2(zeno::format("faceset_{:04}", i), n);
            IFaceSet faceSet = mesh.getFaceSet(n);
            IFaceSetSchema::Sample faceSetSample = faceSet.getSchema().getValue();
            size_t s = faceSetSample.getFaces()->size();
            for (auto j = 0; j < s; j++) {
                int f = faceSetSample.getFaces()->get()[j];
                faceset[f] = i;
            }
        }
    }

    return prim;
}

static std::shared_ptr<PrimitiveObject> foundABCSubd(Alembic::AbcGeom::ISubDSchema &subd, int frameid, bool read_done, bool read_face_set) {
    auto prim = std::make_shared<PrimitiveObject>();

    std::shared_ptr<Alembic::AbcCoreAbstract::v12::TimeSampling> time = subd.getTimeSampling();
    float time_per_cycle =  time->getTimeSamplingType().getTimePerCycle();
    double start = time->getStoredTimes().front();
    int start_frame = std::lround(start / time_per_cycle );
    set_time_info(prim->userData(), time->getTimeSamplingType(), start, int(subd.getNumSamples()));

    int sample_index = clamp(frameid - start_frame, 0, (int)subd.getNumSamples() - 1);
    ISampleSelector iSS = Alembic::Abc::v12::ISampleSelector((Alembic::AbcCoreAbstract::index_t)sample_index);
    Alembic::AbcGeom::ISubDSchema::Sample mesamp = subd.getValue(iSS);

    if (auto marr = mesamp.getPositions()) {
        if (!read_done) {
            log_debug("[alembic] totally {} positions", marr->size());
        }
        auto &parr = prim->verts;
        for (size_t i = 0; i < marr->size(); i++) {
            auto const &val = (*marr)[i];
            parr.emplace_back(val[0], val[1], val[2]);
        }
    }

    read_velocity(prim, mesamp.getVelocities(), read_done);

    if (auto marr = mesamp.getFaceIndices()) {
        if (!read_done) {
            log_debug("[alembic] totally {} face indices", marr->size());
        }
        auto &parr = prim->loops;
        for (size_t i = 0; i < marr->size(); i++) {
            int ind = (*marr)[i];
            parr.push_back(ind);
        }
    }

    if (auto marr = mesamp.getFaceCounts()) {
        if (!read_done) {
            log_debug("[alembic] totally {} faces", marr->size());
        }
        auto &loops = prim->loops;
        auto &parr = prim->polys;
        int base = 0;
        for (size_t i = 0; i < marr->size(); i++) {
            int cnt = (*marr)[i];
            parr.emplace_back(base, cnt);
            base += cnt;
        }
    }
    if (auto uv = subd.getUVsParam()) {
        auto uvsamp =
            uv.getIndexedValue(Alembic::Abc::v12::ISampleSelector((Alembic::AbcCoreAbstract::index_t)sample_index));
        int value_size = (int)uvsamp.getVals()->size();
        int index_size = (int)uvsamp.getIndices()->size();
        if (!read_done) {
            log_debug("[alembic] totally {} uv value", value_size);
            log_debug("[alembic] totally {} uv indices", index_size);
            if (prim->loops.size() == index_size) {
                log_debug("[alembic] uv per face");
            } else if (prim->verts.size() == index_size) {
                log_debug("[alembic] uv per vertex");
            } else {
                log_error("[alembic] error uv indices");
            }
        }
        prim->uvs.resize(value_size);
        {
            auto marr = uvsamp.getVals();
            for (size_t i = 0; i < marr->size(); i++) {
                auto const &val = (*marr)[i];
                prim->uvs[i] = {val[0], val[1]};
            }
        }
        if (prim->loops.size() == index_size) {
            prim->loops.add_attr<int>("uvs");
            for (auto i = 0; i < prim->loops.size(); i++) {
                prim->loops.attr<int>("uvs")[i] = (*uvsamp.getIndices())[i];
            }
        }
        else if (prim->verts.size() == index_size) {
            prim->loops.add_attr<int>("uvs");
            for (auto i = 0; i < prim->loops.size(); i++) {
                prim->loops.attr<int>("uvs")[i] = prim->loops[i];
            }
        }
    }
    if (!prim->loops.has_attr("uvs")) {
        if (!read_done) {
            log_warn("[alembic] Not found uv, auto fill zero.");
        }
        prim->uvs.resize(1);
        prim->uvs[0] = zeno::vec2f(0, 0);
        prim->loops.add_attr<int>("uvs");
        for (auto i = 0; i < prim->loops.size(); i++) {
            prim->loops.attr<int>("uvs")[i] = 0;
        }
    }
    ICompoundProperty arbattrs = subd.getArbGeomParams();
    read_attributes(prim, arbattrs, iSS, read_done);
    ICompoundProperty usrData = subd.getUserProperties();
    read_user_data(prim, usrData, iSS, read_done);

    if (read_face_set) {
        auto &faceset = prim->polys.add_attr<int>("faceset");
        std::fill(faceset.begin(), faceset.end(), -1);
        auto &ud = prim->userData();
        std::vector<std::string> faceSetNames;
        subd.getFaceSetNames(faceSetNames);
        ud.set2("faceset_count", int(faceSetNames.size()));
        for (auto i = 0; i < faceSetNames.size(); i++) {
            auto n = faceSetNames[i];
            ud.set2(zeno::format("faceset_{:04}", i), n);
            IFaceSet faceSet = subd.getFaceSet(n);
            IFaceSetSchema::Sample faceSetSample = faceSet.getSchema().getValue();
            size_t s = faceSetSample.getFaces()->size();
            for (auto j = 0; j < s; j++) {
                int f = faceSetSample.getFaces()->get()[j];
                faceset[f] = i;
            }
        }
    }

    return prim;
}

static std::shared_ptr<CameraInfo> foundABCCamera(Alembic::AbcGeom::ICameraSchema &cam, int frameid) {
    CameraInfo cam_info;
    std::shared_ptr<Alembic::AbcCoreAbstract::v12::TimeSampling> time = cam.getTimeSampling();
    float time_per_cycle =  time->getTimeSamplingType().getTimePerCycle();
    double start = time->getStoredTimes().front();
    int start_frame = std::lround(start / time_per_cycle );
    int sample_index = clamp(frameid - start_frame, 0, (int)cam.getNumSamples() - 1);

    auto samp = cam.getValue(Alembic::Abc::v12::ISampleSelector((Alembic::AbcCoreAbstract::index_t)sample_index));
    cam_info.focal_length = samp.getFocalLength();
    cam_info._near = samp.getNearClippingPlane();
    cam_info._far = samp.getFarClippingPlane();
    cam_info.horizontalAperture = samp.getHorizontalAperture() * 10;
    cam_info.verticalAperture = samp.getVerticalAperture() * 10;
    log_info(
        "[alembic] Camera focal_length: {}, near: {}, far: {}",
        cam_info.focal_length,
        cam_info._near,
        cam_info._far
    );
    return std::make_shared<CameraInfo>(cam_info);
}

static Alembic::Abc::v12::M44d foundABCXform(Alembic::AbcGeom::IXformSchema &xfm, int frameid) {
    std::shared_ptr<Alembic::AbcCoreAbstract::v12::TimeSampling> time = xfm.getTimeSampling();
    float time_per_cycle =  time->getTimeSamplingType().getTimePerCycle();
    double start = time->getStoredTimes().front();
    int start_frame = std::lround(start / time_per_cycle );
    int sample_index = clamp(frameid - start_frame, 0, (int)xfm.getNumSamples() - 1);

    auto samp = xfm.getValue(Alembic::Abc::v12::ISampleSelector((Alembic::AbcCoreAbstract::index_t)sample_index));
    return samp.getMatrix();
}

static std::shared_ptr<PrimitiveObject> foundABCPoints(Alembic::AbcGeom::IPointsSchema &mesh, int frameid, bool read_done) {
    auto prim = std::make_shared<PrimitiveObject>();

    std::shared_ptr<Alembic::AbcCoreAbstract::v12::TimeSampling> time = mesh.getTimeSampling();
    float time_per_cycle =  time->getTimeSamplingType().getTimePerCycle();
    double start = time->getStoredTimes().front();
    int start_frame = std::lround(start / time_per_cycle );
    set_time_info(prim->userData(), time->getTimeSamplingType(), start, int(mesh.getNumSamples()));

    int sample_index = clamp(frameid - start_frame, 0, (int)mesh.getNumSamples() - 1);
    auto iSS = Alembic::Abc::v12::ISampleSelector((Alembic::AbcCoreAbstract::index_t)sample_index);
    Alembic::AbcGeom::IPointsSchema::Sample mesamp = mesh.getValue(iSS);
    if (auto marr = mesamp.getPositions()) {
        if (!read_done) {
            log_info("[alembic] totally {} positions", marr->size());
        }
        auto &parr = prim->verts;
        for (size_t i = 0; i < marr->size(); i++) {
            auto const &val = (*marr)[i];
            parr.emplace_back(val[0], val[1], val[2]);
        }
    }

    {
        auto &ids = prim->verts.add_attr<int>("id");
        auto count = mesamp.getIds()->size();
        for (auto i = 0; i < count; i++) {
            ids[i] = mesamp.getIds()->operator[](i);
        }
    }
    read_velocity(prim, mesamp.getVelocities(), read_done);
    ICompoundProperty arbattrs = mesh.getArbGeomParams();
    read_attributes(prim, arbattrs, iSS, read_done);
    ICompoundProperty usrData = mesh.getUserProperties();
    read_user_data(prim, usrData, iSS, read_done);
    return prim;
}

static std::shared_ptr<PrimitiveObject> foundABCCurves(Alembic::AbcGeom::ICurvesSchema &mesh, int frameid, bool read_done) {
    auto prim = std::make_shared<PrimitiveObject>();

    std::shared_ptr<Alembic::AbcCoreAbstract::v12::TimeSampling> time = mesh.getTimeSampling();
    float time_per_cycle =  time->getTimeSamplingType().getTimePerCycle();
    double start = time->getStoredTimes().front();
    int start_frame = std::lround(start / time_per_cycle );
    set_time_info(prim->userData(), time->getTimeSamplingType(), start, int(mesh.getNumSamples()));

    int sample_index = clamp(frameid - start_frame, 0, (int)mesh.getNumSamples() - 1);
    auto iSS = Alembic::Abc::v12::ISampleSelector((Alembic::AbcCoreAbstract::index_t)sample_index);
    Alembic::AbcGeom::ICurvesSchema::Sample mesamp = mesh.getValue(iSS);
    if (auto marr = mesamp.getPositions()) {
        if (!read_done) {
            log_info("[alembic] totally {} positions", marr->size());
        }
        auto &parr = prim->verts;
        for (size_t i = 0; i < marr->size(); i++) {
            auto const &val = (*marr)[i];
            parr.emplace_back(val[0], val[1], val[2]);
        }
    }
    read_velocity(prim, mesamp.getVelocities(), read_done);
    {
        auto &parr = prim->lines;
        auto numCurves = mesamp.getCurvesNumVertices()->size();
        std::size_t offset = 0;
        for (auto i = 0; i < numCurves; i++) {
            auto count = mesamp.getCurvesNumVertices()->operator[](i);
            for (auto j = 0; j < count-1; j++) {
                parr.push_back(vec2i(offset + j, offset + j + 1));
            }
            offset += count;
        }
    }
    ICompoundProperty arbattrs = mesh.getArbGeomParams();
    read_attributes(prim, arbattrs, iSS, read_done);
    ICompoundProperty usrData = mesh.getUserProperties();
    read_user_data(prim, usrData, iSS, read_done);
    return prim;
}

void traverseABC(
    Alembic::AbcGeom::IObject &obj,
    ABCTree &tree,
    int frameid,
    bool read_done,
    bool read_face_set,
    std::string path
) {
    {
        auto const &md = obj.getMetaData();
        if (!read_done) {
            log_debug("[alembic] meta data: [{}]", md.serialize());
        }
        tree.name = obj.getName();
        path = zeno::format("{}/{}", path, tree.name);

        if (Alembic::AbcGeom::IPolyMesh::matches(md)) {
            if (!read_done) {
                log_debug("[alembic] found a mesh [{}]", obj.getName());
            }

            Alembic::AbcGeom::IPolyMesh meshy(obj);
            auto &mesh = meshy.getSchema();
            tree.prim = foundABCMesh(mesh, frameid, read_done, read_face_set);
            tree.prim->userData().set2("_abc_name", obj.getName());
            tree.prim->userData().set2("_abc_path", path);
        } else if (Alembic::AbcGeom::IXformSchema::matches(md)) {
            if (!read_done) {
                log_debug("[alembic] found a Xform [{}]", obj.getName());
            }
            Alembic::AbcGeom::IXform xfm(obj);
            auto &cam_sch = xfm.getSchema();
            tree.xform = foundABCXform(cam_sch, frameid);
        } else if (Alembic::AbcGeom::ICameraSchema::matches(md)) {
            if (!read_done) {
                log_debug("[alembic] found a Camera [{}]", obj.getName());
            }
            Alembic::AbcGeom::ICamera cam(obj);
            auto &cam_sch = cam.getSchema();
            tree.camera_info = foundABCCamera(cam_sch, frameid);
        } else if(Alembic::AbcGeom::IPointsSchema::matches(md)) {
            if (!read_done) {
                log_debug("[alembic] found points [{}]", obj.getName());
            }
            Alembic::AbcGeom::IPoints points(obj);
            auto &points_sch = points.getSchema();
            tree.prim = foundABCPoints(points_sch, frameid, read_done);
            tree.prim->userData().set2("_abc_name", obj.getName());
            tree.prim->userData().set2("_abc_path", path);
            tree.prim->userData().set2("faceset_count", 0);
        } else if(Alembic::AbcGeom::ICurvesSchema::matches(md)) {
            if (!read_done) {
                log_debug("[alembic] found curves [{}]", obj.getName());
            }
            Alembic::AbcGeom::ICurves curves(obj);
            auto &curves_sch = curves.getSchema();
            tree.prim = foundABCCurves(curves_sch, frameid, read_done);
            tree.prim->userData().set2("_abc_name", obj.getName());
            tree.prim->userData().set2("_abc_path", path);
            tree.prim->userData().set2("faceset_count", 0);
        } else if (Alembic::AbcGeom::ISubDSchema::matches(md)) {
            if (!read_done) {
                log_debug("[alembic] found SubD [{}]", obj.getName());
            }
            Alembic::AbcGeom::ISubD subd(obj);
            auto &subd_sch = subd.getSchema();
            tree.prim = foundABCSubd(subd_sch, frameid, read_done, read_face_set);
            tree.prim->userData().set2("_abc_name", obj.getName());
            tree.prim->userData().set2("_abc_path", path);
        }
    }

    size_t nch = obj.getNumChildren();
    if (!read_done) {
        log_debug("[alembic] found {} children", nch);
    }

    for (size_t i = 0; i < nch; i++) {
        auto const &name = obj.getChildHeader(i).getName();
        if (!read_done) {
            log_debug("[alembic] at {} name: [{}]", i, name);
        }

        Alembic::AbcGeom::IObject child(obj, name);

        auto childTree = std::make_shared<ABCTree>();
        traverseABC(child, *childTree, frameid, read_done, read_face_set, path);
        tree.children.push_back(std::move(childTree));
    }
}

Alembic::AbcGeom::IArchive readABC(std::string const &path) {
    std::string native_path = std::filesystem::u8path(path).string();
    std::string hdr;
    {
        char buf[5];
        std::memset(buf, 0, 5);
        auto fp = std::fopen(native_path.c_str(), "rb");
        if (!fp)
            throw Exception("[alembic] cannot open file for read: " + path);
        std::fread(buf, 4, 1, fp);
        std::fclose(fp);
        hdr = buf;
    }
    if (hdr == "\x89HDF") {
        log_info("[alembic] opening as HDF5 format");
        return {Alembic::AbcCoreHDF5::ReadArchive(), native_path};
    } else if (hdr == "Ogaw") {
        log_info("[alembic] opening as Ogawa format");
        return {Alembic::AbcCoreOgawa::ReadArchive(), native_path};
    } else {
        throw Exception("[alembic] unrecognized ABC header: [" + hdr + "]");
    }
}

struct ReadAlembic : INode {
    Alembic::Abc::v12::IArchive archive;
    std::string usedPath;
    bool read_done = false;
    virtual void apply() override {
        int frameid;
        if (has_input("frameid")) {
            frameid = get_input<NumericObject>("frameid")->get<int>();
        } else {
            frameid = getGlobalState()->frameid;
        }
        auto abctree = std::make_shared<ABCTree>();
        {
            auto path = get_input<StringObject>("path")->get();
            if (usedPath != path) {
                read_done = false;
            }
            if (read_done == false) {
                archive = readABC(path);
            }
            double start, _end;
            GetArchiveStartAndEndTime(archive, start, _end);
            // fmt::print("GetArchiveStartAndEndTime: {}\n", start);
            // fmt::print("archive.getNumTimeSamplings: {}\n", archive.getNumTimeSamplings());
            auto obj = archive.getTop();
            bool read_face_set = get_input2<bool>("read_face_set");
            traverseABC(obj, *abctree, frameid, read_done, read_face_set, "");
            read_done = true;
            usedPath = path;
        }
        {
            auto namelist = std::make_shared<zeno::ListObject>();
            abctree->visitPrims([&] (auto const &p) {
                auto &ud = p->userData();
                auto _abc_path = ud.template get2<std::string>("_abc_path", "");
                namelist->arr.push_back(std::make_shared<StringObject>(_abc_path));
            });
            auto &ud = abctree->userData();
            ud.set2("prim_count", int(namelist->arr.size()));
            for (auto i = 0; i < namelist->arr.size(); i++) {
                auto n = namelist->arr[i];
                ud.set2(zeno::format("path_{:04}", i), n);
            }
            set_output("namelist", namelist);
        }
        set_output("abctree", std::move(abctree));
    }
};

ZENDEFNODE(ReadAlembic, {
    {
        {"readpath", "path"},
        {"bool", "read_face_set", "0"},
        {"frameid"},
    },
    {
        {"ABCTree", "abctree"},
        "namelist",
    },
    {},
    {"alembic"},
});

struct AlembicSplitByName: INode {
    void apply() override {
        auto prim = get_input<PrimitiveObject>("prim");
        int faceset_count = prim->userData().get2<int>("faceset_count");
        {
            auto namelist = std::make_shared<zeno::ListObject>();
            for (auto f = 0; f < faceset_count; f++) {
                auto name = prim->userData().get2<std::string>(zeno::format("faceset_{:04}", f));
                namelist->arr.push_back(std::make_shared<StringObject>(name));
            }
            set_output("namelist", namelist);
        }

        auto dict = std::make_shared<zeno::DictObject>();
        if (prim->polys.size()) {
            std::map<int, std::vector<int>> faceset_map;
            for (auto f = 0; f < faceset_count; f++) {
                faceset_map[f] = {};
            }
            auto &faceset = prim->polys.add_attr<int>("faceset");
            for (auto j = 0; j < faceset.size(); j++) {
                auto f = faceset[j];
                faceset_map[f].push_back(j);
            }
            for (auto f = 0; f < faceset_count; f++) {
                auto name = prim->userData().get2<std::string>(zeno::format("faceset_{:04}", f));
                auto new_prim = std::dynamic_pointer_cast<PrimitiveObject>(prim->clone());
                new_prim->polys.resize(faceset_map[f].size());
                for (auto i = 0; i < faceset_map[f].size(); i++) {
                    new_prim->polys[i] = prim->polys[faceset_map[f][i]];
                }
                new_prim->polys.foreach_attr<AttrAcceptAll>([&](auto const &key, auto &arr) {
                    using T = std::decay_t<decltype(arr[0])>;
                    auto &attr = prim->polys.attr<T>(key);
                    for (auto i = 0; i < arr.size(); i++) {
                        arr[i] = attr[faceset_map[f][i]];
                    }
                });
                new_prim->userData().del("faceset_count");
                for (auto j = 0; j < faceset_count; j++) {
                    new_prim->userData().del(zeno::format("faceset_{:04}", j));
                }
                new_prim->userData().set2("_abc_faceset", name);
                dict->lut[name] = std::move(new_prim);
            }
        }
        else if (prim->tris.size()) {
            std::map<int, std::vector<int>> faceset_map;
            for (auto f = 0; f < faceset_count; f++) {
                faceset_map[f] = {};
            }
            auto &faceset = prim->tris.add_attr<int>("faceset");
            for (auto j = 0; j < faceset.size(); j++) {
                auto f = faceset[j];
                faceset_map[f].push_back(j);
            }
            for (auto f = 0; f < faceset_count; f++) {
                auto name = prim->userData().get2<std::string>(zeno::format("faceset_{:04}", f));
                auto new_prim = std::dynamic_pointer_cast<PrimitiveObject>(prim->clone());
                new_prim->tris.resize(faceset_map[f].size());
                for (auto i = 0; i < faceset_map[f].size(); i++) {
                    new_prim->tris[i] = prim->tris[faceset_map[f][i]];
                }
                new_prim->tris.foreach_attr<AttrAcceptAll>([&](auto const &key, auto &arr) {
                    using T = std::decay_t<decltype(arr[0])>;
                    auto &attr = prim->tris.attr<T>(key);
                    for (auto i = 0; i < arr.size(); i++) {
                        arr[i] = attr[faceset_map[f][i]];
                    }
                });
                new_prim->userData().del("faceset_count");
                for (auto j = 0; j < faceset_count; j++) {
                    new_prim->userData().del(zeno::format("faceset_{:04}", j));
                }
                new_prim->userData().set2("_abc_faceset", name);
                dict->lut[name] = std::move(new_prim);
            }
        }
        set_output("dict", dict);
    }
};

ZENDEFNODE(AlembicSplitByName, {
    {
        {"prim"},
    },
    {
        {"DictObject", "dict"},
        {"ListObject", "namelist"},
    },
    {},
    {"alembic"},
});

struct CopyPosAndNrmByIndex: INode {
    void apply() override {
        auto prim = get_input<PrimitiveObject>("prim");
        auto prims = get_input<ListObject>("list")->get<PrimitiveObject>();
        for (auto p: prims) {
            size_t size = p->size();
            auto index = p->attr<int>("index");
            for (auto i = 0; i < size; i++) {
                prim->verts[index[i]] = p->verts[i];
            }
            if (prim->verts.attr_is<vec3f>("nrm")) {
                auto &nrm = prim->verts.attr<vec3f>("nrm");
                auto &nrm_sub = p->verts.attr<vec3f>("nrm");
                for (auto i = 0; i < size; i++) {
                    nrm[index[i]] = nrm_sub[i];
                }
            }
        }

        set_output("out", prim);
    }
};

ZENDEFNODE(CopyPosAndNrmByIndex, {
    {
        {"prim"},
        {"list", "list"},
    },
    {
        {"out"},
    },
    {},
    {"alembic"},
});

struct PrimsFilterInUserdata: INode {
    void apply() override {
        auto prims = get_input<ListObject>("list")->get<PrimitiveObject>();
        auto filter_str = get_input2<std::string>("filters");
        std::vector<std::string> filters = zeno::split_str(filter_str);
        std::vector<std::string> filters_;
        auto out_list = std::make_shared<ListObject>();

        for (auto &s: filters) {
            if (s.length() > 0) {
                filters_.push_back(s);
            }
        }

        auto name = get_input2<std::string>("name");
        auto contain = get_input2<bool>("contain");
        for (auto p: prims) {
            auto &ud = p->userData();
            bool this_contain = false;
            if (ud.has<std::string>(name)) {
                this_contain = std::count(filters_.begin(), filters_.end(), ud.get2<std::string>(name)) > 0;
            }
            else if (ud.has<int>(name)) {
                this_contain = std::count(filters_.begin(), filters_.end(), std::to_string(ud.get2<int>(name))) > 0;
            }
            else if (ud.has<float>(name)) {
                this_contain = std::count(filters_.begin(), filters_.end(), std::to_string(ud.get2<float>(name))) > 0;
            }
            bool insert = (contain && this_contain) || (!contain && !this_contain);
            if (insert) {
                out_list->arr.push_back(p);
            }
        }
        set_output("out", out_list);
    }
};

ZENDEFNODE(PrimsFilterInUserdata, {
    {
        {"list", "list"},
        {"string", "name", ""},
        {"string", "filters"},
        {"bool", "contain", "1"},
    },
    {
        {"out"},
    },
    {},
    {"alembic"},
});

#ifdef ZENO_WITH_PYTHON3
static PyObject * pycheck(PyObject *pResult) {
    if (pResult == nullptr) {
        PyErr_Print();
        throw zeno::makeError("python err");
    }
    return pResult;
}

static void pycheck(int result) {
    if (result != 0) {
        PyErr_Print();
        throw zeno::makeError("python err");
    }
}
struct PrimsFilterInUserdataPython: INode {
    void apply() override {
        auto prims = get_input<ListObject>("list")->get<PrimitiveObject>();
        auto py_code = get_input2<std::string>("py_code");
        Py_Initialize();
        zeno::scope_exit init_defer([=]{ Py_Finalize(); });

        auto out_list = std::make_shared<ListObject>();
        for (auto p: prims) {
            PyObject* userGlobals = PyDict_New();
            zeno::scope_exit userGlobals_defer([=]{ Py_DECREF(userGlobals); });

            PyObject* innerDict = PyDict_New();
            zeno::scope_exit innerDict_defer([=]{ Py_DECREF(innerDict); });

            auto &ud = p->userData();
            for (auto i = ud.begin(); i != ud.end(); i++) {
                auto key = i->first;
                if (ud.has<std::string>(key)) {
                    auto value = ud.get2<std::string>(key);
                    PyObject* pyInnerValue = PyUnicode_DecodeUTF8(key.c_str(), key.size(), "strict");
                    pycheck(PyDict_SetItemString(innerDict, key.c_str(), pyInnerValue));
                }
                else if (ud.has<float>(key)) {
                    auto value = ud.get2<float>(key);
                    PyObject* pyInnerValue = PyFloat_FromDouble(value);
                    pycheck(PyDict_SetItemString(innerDict, key.c_str(), pyInnerValue));
                }
                else if (ud.has<int>(key)) {
                    auto value = ud.get2<int>(key);
                    PyObject* pyInnerValue = PyLong_FromLong(value);
                    pycheck(PyDict_SetItemString(innerDict, key.c_str(), pyInnerValue));
                }
            }

            PyDict_SetItemString(userGlobals, "ud", innerDict);

            PyObject* pResult = pycheck(PyRun_String(py_code.c_str(), Py_file_input, userGlobals, nullptr));
            zeno::scope_exit pResult_defer([=]{ Py_DECREF(pResult); });
            PyObject* pValue = pycheck(PyRun_String("result", Py_eval_input, userGlobals, nullptr));
            zeno::scope_exit pValue_defer([=]{ Py_DECREF(pValue); });
            int need_insert = PyLong_AsLong(pValue);

            if (need_insert > 0) {
                out_list->arr.push_back(p);
            }
        }
        set_output("out", out_list);
    }
};

ZENDEFNODE(PrimsFilterInUserdataPython, {
    {
        {"list", "list"},
        {"multiline_string", "py_code", "result = len(ud['label']) > 2"},
    },
    {
        {"out"},
    },
    {},
    {"alembic"},
});
#endif


} // namespace zeno
