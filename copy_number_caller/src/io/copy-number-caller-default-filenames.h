#pragma once

namespace xoos::cnc {

// likelihood output
constexpr auto* kDefaultSomaticCNCallsetVcfOutput = "somatic_copy_number_callset.vcf.gz";
constexpr auto* kDefaultGermlineCNCallsetVcfOutput = "germline_copy_number_callset.vcf.gz";
constexpr auto* kDefaultSomaticCNCallsetSegOutput = "somatic_copy_number_callset.seg";
constexpr auto* kDefaultGermlineCNCallsetSegOutput = "germline_copy_number_callset.seg";

constexpr auto* kDefaultBafOutput = "bafs.bed";
constexpr auto* kDefaultBafBwOutput = "bafs.bw";
constexpr auto* kDefaultBafSegOutput = "baf_segments.seg";
constexpr auto* kDefaultLogRsOutput = "log_ratios.bed";
constexpr auto* kDefaultLogRsBwOutput = "log_ratios.bw";
constexpr auto* kDefaultLogRsSegOutput = "log_ratio_segments.seg";

// coverage output
constexpr auto* kDefaultCoverageOutput = "coverage.bed";
constexpr auto* kDefaultCorrectedCoverageOutput = "corrected_coverage.bed";
constexpr auto* kDefaultNormalCoverageOutput = "normal_coverage.bed";
constexpr auto* kDefaultNormalCorrectedCoverageOutput = "normal_corrected_coverage.bed";
constexpr auto* kDefaultTumorCoverageOutput = "tumor_coverage.bed";
constexpr auto* kDefaultTumorCorrectedCoverageOutput = "tumor_corrected_coverage.bed";

// augmented baits output
constexpr auto* kDefaultAugmentedBaitsOutput = "augmented_intervals.bed";

constexpr auto* kDefaultPurityPloidyGridOutput = "purity_ploidy_grid.tsv";
constexpr auto* kDefaultPurityPloidyOutput = "purity_ploidy.tsv";
constexpr auto* kDefaultMapQOutput = "mapping_quality.bed";
constexpr auto* kDefaultMetricsOutput = "metrics.tsv";
constexpr auto* kDefaultIgvXmlOutput = "igv_visualization.xml";
constexpr auto* kDefaultTaskflowGraphOutput = "taskflow.dot";

// denoise output
constexpr auto* kDefaultTumorSexOutput = "tumor_sex_prediction.tsv";
constexpr auto* kDefaultDenoisedLogRsOutput = "denoised_log_ratios.bed";

// panel of normals output
constexpr auto* kDefaultPanelOfNormalsOutput = "panel_of_normals.tsv";
constexpr auto* kDefaultPonAutosomalIntervalMediansOutput = "pon_autosomal_interval_medians.tsv";

// seg to vcf output
constexpr auto* kDefaultVcfOutput = "segments.vcf";

// segmentation output
constexpr auto* kDefaultSegmentationSegOutput = "segments.seg";

}  // namespace xoos::cnc
