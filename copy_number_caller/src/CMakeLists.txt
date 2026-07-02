find_package(Armadillo REQUIRED)
find_package(BigWig REQUIRED)
find_package(CLI11 REQUIRED)
find_package(Taskflow REQUIRED)
find_package(ZLIB REQUIRED)
find_package(xoos-cgranges REQUIRED)
find_package(xoos-hclust-cpp REQUIRED)
find_package(xoos-htslib REQUIRED)
find_package(vincentlaucsb-csv-parser REQUIRED)

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
        baits.cpp
        bigwig_cpp.cpp
        coverage.cpp
        io/fai.cpp
        filter-observations.cpp
        intervals/interval-tools.cpp
        mapq-utils.cpp
        io/write-baf-values.cpp
        io/write-bigwig.cpp
        observations.cpp
        panel-of-normals.cpp
        sex.cpp
        singular-value-decomposition.cpp
        utility/utility-functions.cpp
        io/write-igv-xml.cpp
)

target_include_directories(copy_number_caller_core SYSTEM INTERFACE .)
target_include_directories(copy_number_caller_core PRIVATE .)

target_link_libraries(copy_number_caller_core
        PUBLIC
        armadillo::armadillo
        BigWig::BigWig
        xoos-cgranges::xoos-cgranges
        vincentlaucsb-csv-parser::vincentlaucsb-csv-parser
        PRIVATE
        CLI11::CLI11
        Taskflow::Taskflow
        xoos-hclust-cpp::xoos-hclust-cpp
)

target_link_xoos_libraries(copy_number_caller_core
        PRIVATE
        log
        io
        util
        enum
        error
        gc_correct
        types
        stats
)

target_link_xoos_libraries(copy_number_caller_core
        PUBLIC
        cli
)


target_compile_definitions(copy_number_caller_core
        PRIVATE
        PROGRAM_NAME="copy_number_caller"
        VERSION="${VERSION}"
)

add_default_library(generate_panel_of_normals_core
        generate-panel-of-normals/define-generate-panel-of-normals-options.cpp
        generate-panel-of-normals/generate-panel-of-normals.cpp
        generate-panel-of-normals/generate-panel-of-normals-cli.cpp
)

target_include_directories(generate_panel_of_normals_core SYSTEM INTERFACE .)
target_include_directories(generate_panel_of_normals_core PRIVATE .)
target_link_libraries(generate_panel_of_normals_core
        PUBLIC
        armadillo::armadillo
        PRIVATE
        copy_number_caller_core
        BigWig::BigWig
        CLI11::CLI11
        Taskflow::Taskflow
)

target_link_xoos_libraries(generate_panel_of_normals_core
        PRIVATE
        log
        io
        util
        enum
        gc_correct
)

target_compile_definitions(generate_panel_of_normals_core
        PRIVATE
        PROGRAM_NAME="generate_panel_of_normals"
        VERSION="${VERSION}"
)
