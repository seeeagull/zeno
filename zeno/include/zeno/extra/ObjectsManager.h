#pragma once

#include <zeno/utils/MapStablizer.h>
#include <zeno/utils/PolymorphicMap.h>
#include <zeno/utils/disable_copy.h>
#include <zeno/core/IObject.h>
#include <string>
#include <memory>
#include <map>
#include <set>
#include <optional>

namespace zeno {

struct ObjectsManager {
    zeno::MapStablizer<zeno::PolymorphicMap<std::map<
        std::string, std::shared_ptr<zeno::IObject>>>> objects;

    std::map<std::string, std::shared_ptr<zeno::IObject>> lightObjects;
    bool needUpdateLight = true;

    template <class T = void>
    auto pairs() const {
        return objects.pairs<T>();
    }

    template <class T = void>
    auto pairsShared() const {
        return objects.pairsShared<T>();
    }

    ObjectsManager();
    ~ObjectsManager();
    bool load_objects(std::map<std::string, std::shared_ptr<zeno::IObject>> const &objs);

    ZENO_API void clear_objects();
    ZENO_API std::optional<zeno::IObject*> get(std::string nid);

    //---determine update type accord to objs changes---
    enum RenderType
    {
        UNDEFINED = 0,
        UPDATE_ALL,
        UPDATE_LIGHT_CAMERA,
        UPDATE_MATERIAL
    };
    void determineRenderType(std::map<std::string, std::shared_ptr<zeno::IObject>> const& objs);
    RenderType renderType = UNDEFINED;
    std::map<std::string, int> lastToViewNodesType;
};

}