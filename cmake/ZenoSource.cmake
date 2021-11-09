zeno_glob_recurse(source src *.h *.cpp)
target_sources(zeno PRIVATE ${source})
target_include_directories(zeno PUBLIC src)

if (ZENO_WITH_SYCL)
    message("-- Building Zeno with hipSYCL targets: [${HIPSYCL_TARGETS}]")
    find_package(hipSYCL CONFIG REQUIRED)
    target_compile_definitions(zeno PUBLIC -DZENO_WITH_SYCL)
endif()

if (ZENO_WITH_LEGACY)
    message("-- Building Zeno with Legacy Nodes")
    zeno_glob_recurse(source legacy *.h *.cpp)
    target_include_directories(zeno PUBLIC legacy)
    target_sources(zeno PRIVATE ${source})
endif()

if (ZENO_WITH_BACKWARD)
    message("-- Building Zeno with Stack Traceback")
    target_sources(zeno PRIVATE ${BACKWARD_ENABLE})
    add_backward(zeno)
endif()
