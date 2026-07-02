#pragma once

namespace xoos::cnc::cli_opt_name {

// Common
constexpr auto* kOutputDir = "--output-dir";

// Shared I/O
constexpr auto* kBam = "--bam";
constexpr auto* kTumorBam = "--tumor-bam";
constexpr auto* kNormalBam = "--normal-bam";
constexpr auto* kPanelOfNormals = "--panel-of-normals";
constexpr auto* kIntervals = "--intervals";
constexpr auto* kReference = "--reference";
constexpr auto* kLogRatios = "--log-ratios";
constexpr auto* kSegments = "--segments";
constexpr auto* kBafs = "--bafs";
constexpr auto* kMapqs = "--mapqs";

// AugmentBaits Inputs
constexpr auto* kMappabilityBigwig = "--mappability-bigwig";
constexpr auto* kBlocklistBed = "--blocklist-bed";

// AugmentBaits Other
constexpr auto* kIncludeAltContigs = "--include-alt-contigs";

// AugmentBaits Whole Genome
constexpr auto* kWholeGenome = "--whole-genome";
constexpr auto* kIntervalSize = "--interval-size";
constexpr auto* kMinIntervalMappability = "--min-interval-mappability";
constexpr auto* kMinGcContent = "--min-gc-content";
constexpr auto* kMaxGcContent = "--max-gc-content";

// AugmentBaits Baits Options
constexpr auto* kBed = "--bed";
constexpr auto* kNoOffTargets = "--no-off-targets";
constexpr auto* kOnTargetIntervalSize = "--on-target-interval-size";
constexpr auto* kOffTargetMinWidth = "--off-target-min-width";
constexpr auto* kOffTargetTrimLength = "--off-target-trim-length";
constexpr auto* kOffTargetIntervalSize = "--off-target-interval-size";
constexpr auto* kOnTargetMinMappability = "--on-target-min-mappability";
constexpr auto* kOffTargetMinMappability = "--off-target-min-mappability";
constexpr auto* kSexChromosomeTargetMinMappability = "--sex-chromosome-target-min-mappability";

// CalculateCoverage
constexpr auto* kExcludeFlags = "--exclude-flags";
constexpr auto* kSkipDelAndRefBases = "--skip-del-and-ref-bases";

// GCCorrect
constexpr auto* kFirstPassGcSmoothingFraction = "--first-pass-gc-smoothing-fraction";

// Denoise
constexpr auto* kDenoiseMinTargetLength = "--denoise-min-target-length";
constexpr auto* kDenoiseMinPanelMedianCoverage = "--denoise-min-panel-median-coverage";
constexpr auto* kDenoiseMinPanelMedianAndTumorCoverage = "--denoise-min-panel-median-and-tumor-coverage";
constexpr auto* kDenoiseMinOffTargetFilterFrac = "--denoise-min-off-target-filter-frac";
constexpr auto* kDenoiseDisableFilter = "--denoise-disable-filter";
constexpr auto* kTumorCoverage = "--tumor-coverage";
constexpr auto* kPanelOfNormalsCoveragesList = "--panel-of-normals-coverages-list";

// TwoSampleLogR / Shared coverage
constexpr auto* kTumorCoverageFile = "--tumor-coverage";
constexpr auto* kNormalCoverageFile = "--normal-coverage";
constexpr auto* kCoverageFile = "--coverage";

// Segmentation
constexpr auto* kSeedSegments = "--seed-segments";
constexpr auto* kCbsMethod = "--cbs-method";
constexpr auto* kSegmentationMode = "--segmentation-mode";
constexpr auto* kTruncateOutliers = "--truncate-outliers";
constexpr auto* kMinObservationsPerSegment = "--min-observations-per-segment";
constexpr auto* kHierarchicalPruningHeight = "--hierarchical-pruning-height";
constexpr auto* kEnableSegmentUndoing = "--enable-segment-undoing";
constexpr auto* kInitialSegmentUndoingFactor = "--initial-segment-undoing-factor";
constexpr auto* kMaxSegmentsForUndoing = "--max-segments-for-undoing";
constexpr auto* kSegmentUndoingIncrement = "--segment-undoing-increment";
constexpr auto* kSkipTTestMerging = "--skip-t-test-merging";
constexpr auto* kSkipHierarchicalPruning = "--skip-hierarchical-pruning";
constexpr auto* kEnableDhSegmentUndoing = "--enable-dh-segment-undoing";
constexpr auto* kInitialDhSegmentUndoingFactor = "--initial-dh-segment-undoing-factor";
constexpr auto* kMinTValueToSkipPermutationTesting = "--min-t-value-to-skip-permutation-testing";

// Likelihood
constexpr auto* kMinMapqForCalls = "--min-mapq-for-calls";
constexpr auto* kMinCnvLengthToFlag = "--min-cnv-length-to-flag";

// MergeSegments
constexpr auto* kMinMapq = "--min-mapq";
constexpr auto* kUseMapqsFromFile = "--use-mapqs-from-file";
constexpr auto* kRecalculatePerSegmentData = "--recalculate-per-segment-data";

// VCF parsing
constexpr auto* kVcf = "--vcf";
constexpr auto* kForceEnableSomaticVariantParsing = "--force-enable-somatic-variant-parsing";
constexpr auto* kMinVcfNormalDepth = "--min-vcf-normal-depth";
constexpr auto* kMinNormalBaf = "--min-normal-baf";
constexpr auto* kMaxNormalBaf = "--max-normal-baf";
constexpr auto* kMinVcfTumorDepth = "--min-vcf-tumor-depth";
constexpr auto* kMinTumorBaf = "--min-tumor-baf";
constexpr auto* kMaxTumorBaf = "--max-tumor-baf";
constexpr auto* kAlleleDepthTsv = "--allele-depth-tsv";

// Sample metadata
constexpr auto* kSampleId = "--sample-id";
constexpr auto* kSex = "--sex";
constexpr auto* kSampleName = "--sample-name";
constexpr auto* kTumorSampleName = "--tumor-sample-name";
constexpr auto* kNormalSampleName = "--normal-sample-name";
constexpr auto* kTumorPurity = "--tumor-purity";
constexpr auto* kTumorPloidy = "--tumor-ploidy";

// PurityPloidySearch
constexpr auto* kMaxSegmentLogRatioForPurityPloidySearch = "--max-segment-log-ratio-for-purity-ploidy-search";
constexpr auto* kMinSnpsPerSegmentForPurityPloidySearch = "--min-snps-per-segment-for-purity-ploidy-search";

// Coverage thresholds
constexpr auto* kMinTumorCoverage = "--min-tumor-coverage";
constexpr auto* kMinNormalCoverage = "--min-normal-coverage";

// New hidden output-control flags
constexpr auto* kSaveTaskflowGraph = "--save-taskflow-graph";
constexpr auto* kSaveLogRatioSegments = "--save-log-ratio-segments";
constexpr auto* kSaveBafSegments = "--save-baf-segments";

// Performance
constexpr auto* kThreads = "--threads";

// CLI option groups
constexpr auto* kInputOptionsGroup = "Input Options";
constexpr auto* kOutputOptionsGroup = "Output Options";
constexpr auto* kReadFilteringOptionsGroup = "Read Filtering Options";
constexpr auto* kSegmentationOptionsGroup = "Segmentation Options";
constexpr auto* kRegionOptionsGroup = "Region Options";
constexpr auto* kSampleMetadataOptionsGroup = "Sample Metadata Options";
constexpr auto* kVcfParsingOptionsGroup = "VCF Parsing Options";
constexpr auto* kPerformanceOptionsGroup = "Performance Options";
constexpr auto* kPurityPloidySearchOptionsGroup = "Purity/Ploidy Search Options";
constexpr auto* kDenoiseOptionsGroup = "Denoise Options";
constexpr auto* kGCCorrectOptionsGroup = "GC Correct Options";
constexpr auto* kCalculateCoverageOptionsGroup = "Calculate Coverage Options";
constexpr auto* kHiddenOptionsGroup = "";
constexpr auto* kAdvancedOptionsGroup = "Advanced Options";

}  // namespace xoos::cnc::cli_opt_name
