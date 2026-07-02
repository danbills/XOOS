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

find_package(vincentlaucsb-csv-parser REQUIRED)
target_link_libraries(read_collapser_rescue_secondary_core PRIVATE vincentlaucsb-csv-parser::vincentlaucsb-csv-parser)
target_link_libraries(read_collapser_rescue_secondary_core PRIVATE read_collapser_metrics_format_util_core)

target_include_directories(read_collapser_core SYSTEM INTERFACE .)
target_include_directories(read_collapser_core PRIVATE .)

include_directories(
        SYSTEM
        INTERFACE
        ../external/samtools
)
include_directories(
        PRIVATE
        ../external/samtools
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

find_package(xoos-spoa REQUIRED)
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
