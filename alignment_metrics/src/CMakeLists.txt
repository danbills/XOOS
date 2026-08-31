# ============================================================================
# XOOS Alignment Metrics Core Library CMakeLists.txt
# ============================================================================

find_package(xoos-htslib QUIET)
find_package(vincentlaucsb-csv-parser QUIET)

if (xoos-htslib_FOUND AND vincentlaucsb-csv-parser_FOUND)
    add_default_library(alignment_metrics_core
        core/alignment.cpp
        core/hp-error-for-read.cpp
        core/interval.cpp
        core/pileup.cpp
        core/region-lookup.cpp
        core/substitution-lookup.cpp
        core/super-region.cpp
        io/bed-reader.cpp
        metadata/alignment-metadata.cpp
        metadata/dataset-metadata.cpp
        metrics/accuracy-metrics/base-level-accuracy-summary.cpp
        metrics/accuracy-metrics/errors-by-cluster-size.cpp
        metrics/accuracy-metrics/errors-by-read-type.cpp
        metrics/accuracy-metrics/errors-by-substitution-type.cpp
        metrics/accuracy-metrics/hp-stat.cpp
        metrics/accuracy-metrics/qscore-stats.cpp
        metrics/coverage-metrics/coverage-distribution-summary.cpp
        metrics/coverage-metrics/coverage-histograms.cpp
        metrics/coverage-metrics/coverage-uniformity-metrics.cpp
        metrics/read-metrics/read-counts-by-cluster-size.cpp
        metrics/read-metrics/read-length-histograms.cpp
        metrics/read-metrics/read-metrics-summary.cpp
        metrics/read-metrics/read-metrics-te.cpp
        metrics/metrics.cpp
        util/format-util.cpp
        util/homopolymer-util.cpp
        util/region-util.cpp
        alignment-metrics-cli.cpp
        alignment-metrics-options.cpp
        alignment-metrics.cpp
        metrics-worker.cpp
        run-alignment-metrics.cpp
    )
    target_include_directories(alignment_metrics_core
        PRIVATE
        .
    )
    target_include_directories(alignment_metrics_core
        SYSTEM
        INTERFACE
        .
    )
    target_link_xoos_libraries(alignment_metrics_core
        PUBLIC
        histogram
        io
        types
        yc_decode
        PRIVATE
        cli
        enum
        error
        log
        util
    )

    find_package(Taskflow QUIET)
    find_package(Microsoft.GSL QUIET)
    find_package(unordered_dense QUIET)

    if (Taskflow_FOUND AND Microsoft.GSL_FOUND AND unordered_dense_FOUND)
        target_link_libraries(alignment_metrics_core
            PUBLIC
            Taskflow::Taskflow
            PRIVATE
            vincentlaucsb-csv-parser::vincentlaucsb-csv-parser
            Microsoft.GSL::GSL
            unordered_dense::unordered_dense
        )
    endif ()
endif ()

# CUDA GPU Acceleration for alignment_metrics
include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS "CUDA Acceleration: ENABLED for alignment_metrics_core (Compiler: ${CMAKE_CUDA_COMPILER})")

    add_default_library(alignment_metrics_cuda_core
        cuda/alignment_metrics_cuda.cu
    )
    target_include_directories(alignment_metrics_cuda_core PRIVATE . cuda)
    target_include_directories(alignment_metrics_cuda_core SYSTEM INTERFACE . cuda)
    target_link_libraries(alignment_metrics_cuda_core PUBLIC CUDA::cudart)
    target_link_xoos_libraries(alignment_metrics_cuda_core PUBLIC types)
    target_compile_definitions(alignment_metrics_cuda_core PUBLIC XOOS_CUDA_ENABLED=1)
endif ()
