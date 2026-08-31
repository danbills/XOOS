# ============================================================================
# XOOS Read Collapser Executables CMakeLists.txt
# ============================================================================

if (TARGET read_collapser_core)
    add_default_executable(read_collapser
        read-collapser.cpp
    )
    target_link_libraries(read_collapser
        PRIVATE
        read_collapser_core
    )
    target_link_xoos_libraries(read_collapser
        PRIVATE
        io
        util
        cli
        enum
        log
    )
endif ()

if (TARGET read_collapser_rescue_secondary_core)
    add_default_executable(read_collapser_rescue_secondary
        read-collapser-rescue-secondary.cpp
    )
    target_link_libraries(read_collapser_rescue_secondary
        PRIVATE
        read_collapser_rescue_secondary_core
    )
    target_link_xoos_libraries(read_collapser_rescue_secondary
        PRIVATE
        cli
        io
        log
        util
    )
endif ()

if (CMAKE_CUDA_COMPILER)
    add_default_executable(read_collapser_cuda_bench
        read-collapser-cuda-bench.cpp
    )
    target_link_libraries(read_collapser_cuda_bench
        PRIVATE
        read_collapser_cuda_core
    )
endif ()
