# ============================================================================
# XOOS Read Collapser Core Library CMakeLists.txt
# ============================================================================

find_package(xoos-htslib QUIET)
find_package(xoos-spoa QUIET)

if (xoos-htslib_FOUND AND xoos-spoa_FOUND)
    add_default_library(read_collapser_core
        clustering/cluster-downsample.cpp
        clustering/clustering.cpp
        consensus/base-encoder.cpp
        consensus/consensus-matrix.cpp
        consensus/consensus-matrix-builder.cpp
        consensus/consensus-result.cpp
        consensus/consensus.cpp
        consensus/majority-voting.cpp
        core/region-lookup.cpp
        core/region.cpp
        io/alignment-io.cpp
        io/alignment.cpp
        core/region-partition.cpp
        io/fastq-writer.cpp
        mark-duplicate/mark-duplicate.cpp
        metrics/metrics.cpp
        core/read-collapser-cli.cpp
        util/cli-option-util.cpp
        mark-duplicate/mean-base-quality.cpp
        consensus/poa-msa.cpp
        consensus/qscore-calculator.cpp
        util/duplex-util.cpp
        util/read-util.cpp
        util/softclip-util.cpp
        ../external/samtools/bam_cat.c
        util/read-name-util.cpp
    )

    add_default_library(read_collapser_metrics_format_util_core
        metrics/metrics-format-util.cpp
    )

    target_include_directories(read_collapser_metrics_format_util_core SYSTEM INTERFACE .)
    target_include_directories(read_collapser_metrics_format_util_core PRIVATE .)

    find_package(fmt REQUIRED)
    target_link_libraries(read_collapser_metrics_format_util_core PRIVATE fmt::fmt)
    target_link_xoos_libraries(read_collapser_metrics_format_util_core
        PUBLIC
        types
    )

    add_default_library(read_collapser_rescue_secondary_core
        rescue-secondary/rescue-secondary.cpp
        rescue-secondary/rescue-secondary-cli.cpp
        rescue-secondary/rescue-secondary-metrics.cpp
    )

    target_include_directories(read_collapser_rescue_secondary_core SYSTEM INTERFACE .)
    target_include_directories(read_collapser_rescue_secondary_core PRIVATE .)

    target_link_xoos_libraries(read_collapser_rescue_secondary_core
        PRIVATE
        cli
        error
        io
        log
        util
    )

    target_include_directories(read_collapser_core
        PUBLIC
        ../external/samtools
    )

    target_include_directories(read_collapser_core
        SYSTEM
        INTERFACE
        .
    )

    target_include_directories(read_collapser_core
        PRIVATE
        .
    )

    target_link_xoos_libraries(read_collapser_core
        PUBLIC
        concurrent
        PRIVATE
        cli
        enum
        error
        histogram
        io
        log
        util
        yc_decode
    )

    find_package(Boost REQUIRED)
    target_link_libraries(read_collapser_core PRIVATE Boost::headers)
    target_link_libraries(read_collapser_core PRIVATE xoos-spoa::xoos-spoa)

    find_package(Taskflow REQUIRED)
    target_link_libraries(read_collapser_core PRIVATE Taskflow::Taskflow)

    find_package(Microsoft.GSL REQUIRED)
    target_link_libraries(read_collapser_core PUBLIC Microsoft.GSL::GSL)

    find_package(vincentlaucsb-csv-parser REQUIRED)
    target_link_libraries(read_collapser_core PRIVATE vincentlaucsb-csv-parser::vincentlaucsb-csv-parser)

    find_package(ZLIB REQUIRED)
    target_link_libraries(read_collapser_core PRIVATE ZLIB::ZLIB)
    target_link_libraries(read_collapser_core PRIVATE read_collapser_metrics_format_util_core)
endif ()

# CUDA GPU Acceleration for read_collapser
include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS "CUDA Acceleration: ENABLED for read_collapser_core (Compiler: ${CMAKE_CUDA_COMPILER})")

    add_default_library(read_collapser_cuda_core
        cuda/read_collapser_cuda.cu
    )
    target_include_directories(read_collapser_cuda_core PRIVATE . cuda)
    target_include_directories(read_collapser_cuda_core SYSTEM INTERFACE . cuda)
    target_link_libraries(read_collapser_cuda_core PUBLIC CUDA::cudart)
    target_link_xoos_libraries(read_collapser_cuda_core PUBLIC types)
    target_compile_definitions(read_collapser_cuda_core PUBLIC XOOS_CUDA_ENABLED=1)
endif ()
