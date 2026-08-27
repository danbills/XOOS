# ============================================================================
# XOOS Alignment Metrics Executables CMakeLists.txt
# ============================================================================

find_package(vincentlaucsb-csv-parser QUIET)

if (vincentlaucsb-csv-parser_FOUND AND TARGET alignment_metrics_core)
    add_default_executable(alignment_metrics
        alignment-metrics.cpp
    )
    target_link_xoos_libraries(alignment_metrics
        PRIVATE
        cli
        io
        util
        types
        log
    )
    target_link_libraries(alignment_metrics
        PRIVATE
        alignment_metrics_core
        vincentlaucsb-csv-parser::vincentlaucsb-csv-parser
    )
endif ()

if (CMAKE_CUDA_COMPILER)
    add_default_executable(alignment_metrics_cuda_bench
        alignment-metrics-cuda-bench.cpp
    )
    target_link_libraries(alignment_metrics_cuda_bench
        PRIVATE
        alignment_metrics_cuda_core
    )
endif ()
