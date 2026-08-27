# ============================================================================
# XOOS Copy Number Caller Core Library CMakeLists.txt
# ============================================================================

find_package(Armadillo QUIET)
find_package(BigWig QUIET)
find_package(CLI11 QUIET)
find_package(Taskflow QUIET)
find_package(ZLIB QUIET)
find_package(xoos-cgranges QUIET)
find_package(xoos-hclust-cpp QUIET)
find_package(xoos-htslib QUIET)
find_package(vincentlaucsb-csv-parser QUIET)

if (Armadillo_FOUND AND BigWig_FOUND AND CLI11_FOUND AND Taskflow_FOUND AND xoos-htslib_FOUND)
    add_default_library(copy_number_caller_core
        augment-baits/augment-baits.cpp
        calculate-coverage/calculate-coverage-taskflow-graph.cpp
        calculate-coverage/calculate-coverage.cpp
        denoise/denoise.cpp
        gc-correct/gc-correct.cpp
        likelihood/flag-copy-number-calls.cpp
        likelihood/likelihood-flags.cpp
        likelihood/likelihood.cpp
        likelihood/baf-tools.cpp
        merge-segments/merge-segments.cpp
        merge-segments/merge-segments-main.cpp
        misc/define-thread-options.cpp
        one-sample-logr/one-sample-logr.cpp
        one-sample-logr/self-normalize.cpp
        predict-sex/define-predict-sex-options.cpp
        predict-sex/predict-sex-cli.cpp
        predict-sex/predict-sex.cpp
        purity-ploidy-search/purity-ploidy-search.cpp
        copy-number-caller/define-copy-number-caller-options.cpp
        copy-number-caller/copy-number-caller-cli.cpp
        copy-number-caller/common/bam-to-corrected-coverage.cpp
        copy-number-caller/common/check-chr-y-par.cpp
        copy-number-caller/common/coverage-check.cpp
        copy-number-caller/common/write-empty-outputs.cpp
        io/write-logr.cpp
        io/baits-io.cpp
        copy-number-caller/germline-normal-wgs/germline-normal-wgs.cpp
        copy-number-caller/somatic-tumor-normal-wgs/somatic-tumor-normal-wgs.cpp
        io/write-metrics.cpp
        io/column-names.cpp
        seg-to-vcf/seg-to-vcf.cpp
        segmentation/cbs-util.cpp
        segmentation/cbs.cpp
        segmentation/cbs-wrapper.cpp
        segmentation/cbs-for-somatic.cpp
        segmentation/genomic-segments.cpp
        segmentation/interval-trees.cpp
        segmentation/max-t.cpp
        segmentation/parent-specific-binary-segmentation-taskflow-graph.cpp
        segmentation/pscbs.cpp
        segmentation/segment-by-segments-taskflow-graph.cpp
        segmentation/segment-pruning.cpp
        segmentation/segmentation-mode.cpp
        io/write-segments.cpp
        segmentation/read-segments.cpp
        segmentation/undo-segments.cpp
        two-sample-logr/calculate-logrs.cpp
        two-sample-logr/tumor-normal.cpp
        two-sample-logr/two-sample-logr.cpp
        utility/utility.cpp
        utility/vector-math.cpp
        utility/split-string.cpp
        utility/define-logr-options.cpp
        baits.cpp
        bigwig_cpp.cpp
        coverage.cpp
        filter-observations.cpp
        mapq-utils.cpp
        observations.cpp
        panel-of-normals.cpp
        sex.cpp
        singular-value-decomposition.cpp
    )

    target_include_directories(copy_number_caller_core SYSTEM INTERFACE .)
    target_include_directories(copy_number_caller_core PRIVATE .)

    target_link_libraries(copy_number_caller_core
        PUBLIC
        armadillo::armadillo
        BigWig::BigWig
        CLI11::CLI11
        Taskflow::Taskflow
        ZLIB::ZLIB
        xoos-cgranges::xoos-cgranges
        xoos-hclust-cpp::xoos-hclust-cpp
        xoos-htslib::xoos-htslib
        vincentlaucsb-csv-parser::vincentlaucsb-csv-parser
    )

    target_link_xoos_libraries(copy_number_caller_core
        PRIVATE
        log
        io
        util
        enum
        gc_correct
        types
        stats
        PUBLIC
        cli
    )
endif ()

# CUDA GPU Acceleration for copy_number_caller
include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS "CUDA Acceleration: ENABLED for copy_number_caller_core (Compiler: ${CMAKE_CUDA_COMPILER})")

    add_default_library(copy_number_caller_cuda_core
        cuda/cnv_caller_cuda.cu
    )
    target_include_directories(copy_number_caller_cuda_core PRIVATE . cuda)
    target_include_directories(copy_number_caller_cuda_core SYSTEM INTERFACE . cuda)
    target_link_libraries(copy_number_caller_cuda_core PUBLIC CUDA::cudart)
    target_link_xoos_libraries(copy_number_caller_cuda_core PUBLIC types)
    target_compile_definitions(copy_number_caller_cuda_core PUBLIC XOOS_CUDA_ENABLED=1)
endif ()
