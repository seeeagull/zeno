#include <zeno/zeno.h>
#include <zeno/ParticlesObject.h>
#include <zeno/PrimitiveObject.h>
#include <zeno/VDBGrid.h>
#include "tbb/concurrent_vector.h"
#include "tbb/parallel_for.h"
#include "tbb/scalable_allocator.h"
namespace zeno {

struct GetVDBPoints : zeno::INode {
  virtual void apply() override {
    auto grid = get_input("grid")->as<VDBPointsGrid>()->m_grid;

    std::vector<openvdb::points::PointDataTree::LeafNodeType*> leafs;
    grid->tree().getNodes(leafs);
    printf("GetVDBPoints: particle leaf nodes: %d\n", leafs.size());

    auto transform = grid->transformPtr();

    auto ret = zeno::IObject::make<ParticlesObject>();

    for (auto const &leaf: leafs) {
      //attributes
      // Attribute reader
      // Extract the position attribute from the leaf by name (P is position).
      openvdb::points::AttributeArray& positionArray =
        leaf->attributeArray("P");
      // Extract the velocity attribute from the leaf by name (v is velocity).
      openvdb::points::AttributeArray& velocityArray =
        leaf->attributeArray("v");

      using PositionCodec = openvdb::points::FixedPointCodec</*one byte*/false>;
      using VelocityCodec = openvdb::points::TruncateCodec;
      // Create read handles for position and velocity
      openvdb::points::AttributeHandle<openvdb::Vec3f, PositionCodec> positionHandle(positionArray);
      openvdb::points::AttributeHandle<openvdb::Vec3f, VelocityCodec> velocityHandle(velocityArray);

      for (auto iter = leaf->beginIndexOn(); iter; ++iter) {
        openvdb::Vec3R p = positionHandle.get(*iter);
        p += iter.getCoord().asVec3d();
        // https://people.cs.clemson.edu/~jtessen/cpsc8190/OpenVDB-dpawiki.pdf
        p = transform->indexToWorld(p);
        openvdb::Vec3R v = velocityHandle.get(*iter);
        ret->pos.push_back(glm::vec3(p[0], p[1], p[2]));
        ret->vel.push_back(glm::vec3(v[0], v[1], v[2]));
      }
    }
    set_output("pars", ret);
  }
};

static int defGetVDBPoints = zeno::defNodeClass<GetVDBPoints>("GetVDBPoints",
    { /* inputs: */ {
        "grid",
    }, /* outputs: */ {
        "pars",
    }, /* params: */ {
    }, /* category: */ {
      "openvdb",
    }});


//TODO: parallelize using concurrent vector
struct VDBPointsToPrimitive : zeno::INode {
  virtual void apply() override {
    auto grid = get_input("grid")->as<VDBPointsGrid>()->m_grid;

    std::vector<openvdb::points::PointDataTree::LeafNodeType*> leafs;
    grid->tree().getNodes(leafs);
    printf("GetVDBPoints: particle leaf nodes: %d\n", leafs.size());

    auto transform = grid->transformPtr();

    auto ret = zeno::IObject::make<zeno::PrimitiveObject>();
    auto &retpos = ret->add_attr<zeno::vec3f>("pos");
    auto &retvel = ret->add_attr<zeno::vec3f>("vel");

    tbb::concurrent_vector<std::tuple<zeno::vec3f,zeno::vec3f>> data(0);
    tbb::parallel_for((size_t)0, (size_t)leafs.size(), (size_t)1, [&](size_t index)
    //for (auto const &leaf: leafs)
    {
      auto &leaf = leafs[index];
      //attributes
      // Attribute reader
      // Extract the position attribute from the leaf by name (P is position).
      openvdb::points::AttributeArray& positionArray =
        leaf->attributeArray("P");
      // Extract the velocity attribute from the leaf by name (v is velocity).
      openvdb::points::AttributeArray& velocityArray =
        leaf->attributeArray("v");

      using PositionCodec = openvdb::points::FixedPointCodec</*one byte*/false>;
      using VelocityCodec = openvdb::points::TruncateCodec;
      // Create read handles for position and velocity
      openvdb::points::AttributeHandle<openvdb::Vec3f, PositionCodec> positionHandle(positionArray);
      openvdb::points::AttributeHandle<openvdb::Vec3f, VelocityCodec> velocityHandle(velocityArray);

      for (auto iter = leaf->beginIndexOn(); iter; ++iter) {
        openvdb::Vec3R p = positionHandle.get(*iter);
        p += iter.getCoord().asVec3d();
        // https://people.cs.clemson.edu/~jtessen/cpsc8190/OpenVDB-dpawiki.pdf
        p = transform->indexToWorld(p);
        openvdb::Vec3R v = velocityHandle.get(*iter);
        //retpos.emplace_back(p[0], p[1], p[2]);
        //retvel.emplace_back(v[0], v[1], v[2]);
        data.emplace_back(std::make_tuple(zeno::vec3f(p[0],p[1],p[2]), zeno::vec3f(v[0],v[1],v[2])));
      }
    });
    ret->resize(data.size());
    tbb::parallel_for((size_t)0, (size_t)ret->size(), (size_t)1, 
    [&](size_t index)
    {
      retpos[index] = std::get<0>(data[index]);
      retvel[index] = std::get<1>(data[index]);
    });
    set_output("prim", ret);
  }
};

static int defVDBPointsToPrimitive = zeno::defNodeClass<VDBPointsToPrimitive>("VDBPointsToPrimitive",
    { /* inputs: */ {
        "grid",
    }, /* outputs: */ {
        "prim",
    }, /* params: */ {
    }, /* category: */ {
      "openvdb",
    }});

}
