#include "copy-number-caller/define-copy-number-caller-options.h"

#include <xoos/cli/bam-option-util.h>
#include <xoos/cli/fasta-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/bed-validator.h>
#include <xoos/cli/validators/file-extension-validator.h>
#include <xoos/cli/validators/file-permission-validator.h>
#include <xoos/cli/validators/indexed-vcf-validator.h>

#include "copy-number-caller/copy-number-caller-cli-option-names.h"
#include "segmentation/segmentation-mode.h"

namespace xoos::cnc {

namespace opt = cli_opt_name;

// throwaway for CLI option callback that needs to parse the BAM index path from the BAM path
namespace {
fs::path discarded_bam_index{};
}

// =============================================================================
// Add*Option helper implementations
// =============================================================================

// Shared End-to-End options

CLI::Option* AddPanelOfNormalsOption(const cli::AppPtr& app, const std::string& group_name, fs::path& binding) {
  return app->add_option(opt::kPanelOfNormals, binding, "pre-generated panel of normals in HDF5 format")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddIntervalsOption(const cli::AppPtr& app,
                                const std::string& group_name,
                                std::optional<fs::path>& binding) {
  return app->add_option(opt::kIntervals, binding, "Annotated intervals file.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

// Shared IO options

CLI::Option* AddSegmentsOption(const cli::AppPtr& app,
                               const std::string& group_name,
                               std::optional<fs::path>& binding) {
  return app
      ->add_option(
          opt::kSegments, binding, "Initial copy number calls in SEG format, as output by the copy_number_caller.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddLogrsOption(const cli::AppPtr& app, const std::string& group_name, std::optional<fs::path>& binding) {
  return app
      ->add_option(opt::kLogRatios,
                   binding,
                   "Log ratios file in BED format (Contig, Start, End, LogR), as output by the copy_number_caller.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddBafsOption(const cli::AppPtr& app, const std::string& group_name, std::optional<fs::path>& binding) {
  return app
      ->add_option(opt::kBafs,
                   binding,
                   "Tumor BAF values in BED format (Contig, Start, End, BAF, Ref_AD, Alt_AD), as output by the "
                   "copy_number_caller. Used for estimation of allele-specific copy number.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddMapqsOption(const cli::AppPtr& app, const std::string& group_name, std::optional<fs::path>& binding) {
  return app
      ->add_option(opt::kMapqs,
                   binding,
                   "Mapping quality file in BED format (Contig, Start, End, MeanMapQ). Use with "
                   "--mapq-cutoff-to-flag-calls.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddVcfOption(const cli::AppPtr& app,
                          const std::string& group_name,
                          std::optional<fs::path>& binding,
                          const std::string& description) {
  return app->add_option(opt::kVcf, binding, description)
      ->check(cli::FileReadableValidator())
      ->check(cli::FileExtensionValidator({".vcf", ".vcf.gz"}))
      ->group(group_name);
}

// VCF Parsing

CLI::Option* AddNormalSampleMinDepthOption(const cli::AppPtr& app,
                                           const std::string& group_name,
                                           s32& binding,
                                           const s32 default_depth) {
  return app
      ->add_option(opt::kMinVcfNormalDepth,
                   binding,
                   "Min normal depth for extracting variants from a matched tumor-normal or normal-only VCF. "
                   "Variants with depth less than this value are excluded.")
      ->check(CLI::TypeValidator<s32>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(default_depth)
      ->needs(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddNormalSampleMinBafOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app
      ->add_option(opt::kMinNormalBaf,
                   binding,
                   "Minimum BAF for a variant in the normal sample to be considered heterozygous. Variants below "
                   "this value are excluded.")
      ->check(CLI::TypeValidator<float>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(kVcfParsingOptionsDefaultSomaticNormalSampleMinBAF)
      ->needs(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddNormalSampleMaxBafOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app
      ->add_option(opt::kMaxNormalBaf,
                   binding,
                   "Maximum BAF for a variant in the normal sample to be considered heterozygous. Variants above "
                   "this value are excluded.")
      ->check(CLI::TypeValidator<float>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(kVcfParsingOptionsDefaultSomaticNormalSampleMaxBAF)
      ->needs(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddTumorSampleMinDepthOption(const cli::AppPtr& app, const std::string& group_name, s32& binding) {
  return app->add_option(opt::kMinVcfTumorDepth, binding, "Minimum tumor depth for extracting variants.")
      ->check(CLI::TypeValidator<s32>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(kVcfParsingOptionsDefaultSomaticTumorSampleMinDepth)
      ->needs(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddTumorSampleMinBafOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app
      ->add_option(opt::kMinTumorBaf,
                   binding,
                   "Minimum BAF for a variant in the tumor sample to be extracted. Variants below this value are "
                   "excluded.")
      ->check(CLI::TypeValidator<float>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(kVcfParsingOptionsDefaultSomaticTumorSampleMinBAF)
      ->needs(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddTumorSampleMaxBafOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app
      ->add_option(opt::kMaxTumorBaf,
                   binding,
                   "Maximum BAF for a variant in the tumor sample to be extracted. Variants above this value are "
                   "excluded.")
      ->check(CLI::TypeValidator<float>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(kVcfParsingOptionsDefaultSomaticTumorSampleMaxBAF)
      ->needs(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddAdInOption(const cli::AppPtr& app, const std::string& group_name, std::optional<fs::path>& binding) {
  return app
      ->add_option(
          opt::kAlleleDepthTsv, binding, "Alternative to VCF. Requires Chromosome, Position, and TumorBaf columns.")
      ->excludes(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddForceEnableSomaticVariantParsingOption(const cli::AppPtr& app,
                                                       const std::string& group_name,
                                                       bool& binding) {
  return app
      ->add_flag(opt::kForceEnableSomaticVariantParsing,
                 binding,
                 "Force enable parsing of somatic variants even if the fraction of germline variants is low. Only "
                 "applicable for input VCF that contains flagged somatic variants.")
      ->needs(opt::kVcf)
      ->default_val(false)
      ->group(group_name);
}

// AugmentBaits Inputs
CLI::Option* AddMappabilityBigwigOption(const cli::AppPtr& app, const std::string& group_name, fs::path& binding) {
  return app
      ->add_option(
          opt::kMappabilityBigwig, binding, "Path to an input BigWig file describing mappability across the genome.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddBlocklistBedOption(const cli::AppPtr& app,
                                   const std::string& group_name,
                                   std::optional<fs::path>& binding) {
  return app
      ->add_option(opt::kBlocklistBed,
                   binding,
                   "Path to a BED file of 0-based regions to exclude from the analysis before generating intervals.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

// AugmentBaits Other

CLI::Option* AddIncludeAltContigsOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app
      ->add_flag(opt::kIncludeAltContigs,
                 binding,
                 "include alt contigs when generating the bait file (default: ignore alt and decoy contigs)")
      ->default_val(kAugmentBaitsDefaultIncludeAltContigs)
      ->group(group_name);
}

// AugmentBaits Whole Genome

CLI::Option* AddWholeGenomeOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app
      ->add_flag(opt::kWholeGenome,
                 binding,
                 "generate baits of `--interval-size` across the whole genome. Using this option will "
                 "turn off off-target generation.")
      ->default_val(kAugmentBaitsDefaultWholeGenome)
      ->group(group_name);
}

CLI::Option* AddWholeGenomeIntervalSizeOption(const cli::AppPtr& app,
                                              const std::string& group_name,
                                              size_t& binding,
                                              const size_t default_size) {
  return app
      ->add_option(opt::kIntervalSize,
                   binding,
                   "The size (in base pairs) for each target region of the genome for which log ratios are "
                   "calculated.")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(default_size)
      ->group(group_name);
}

CLI::Option* AddMinGcContentOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app
      ->add_option(opt::kMinGcContent,
                   binding,
                   "Minimum GC content threshold for filtering targets. "
                   "Targets with GC content below this value will be removed.")
      ->default_val(kAugmentBaitsDefaultMinGcContent)
      ->check(CLI::Range(0.0, 1.0))
      ->group(group_name);
}

CLI::Option* AddMaxGcContentOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app
      ->add_option(opt::kMaxGcContent,
                   binding,
                   "Maximum GC content threshold for filtering targets. "
                   "Targets with GC content above this value will be removed.")
      ->default_val(kAugmentBaitsDefaultMaxGcContent)
      ->check(CLI::Range(0.0, 1.0))
      ->group(group_name);
}

CLI::Option* AddWholeGenomeMinMappabilityOption(const cli::AppPtr& app,
                                                const std::string& group_name,
                                                f64& binding,
                                                const f64 default_mappability) {
  return app
      ->add_option(opt::kMinIntervalMappability,
                   binding,
                   "Minimum allowable mappability score for an interval. Intervals with a score less than this "
                   "value are excluded to ensure high-quality read alignment. Mappability scores come from "
                   "`--mappability-bigwig`.")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(default_mappability)
      ->group(group_name);
}

// AugmentBaits Baits Options

CLI::Option* AddBedOption(const cli::AppPtr& app, const std::string& group_name, std::optional<fs::path>& binding) {
  return app->add_option(opt::kBed, binding, "baits file to augment, in BED format")
      ->check(cli::BedFileValidator())
      ->group(group_name);
}

CLI::Option* AddNoOffTargetsOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app->add_flag(opt::kNoOffTargets, binding, "do not generate off-target intervals for baits")
      ->default_val(kAugmentBaitsDefaultNoOffTargets)
      ->group(group_name);
}

CLI::Option* AddOnTargetIntervalSizeOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding) {
  return app
      ->add_option(opt::kOnTargetIntervalSize,
                   binding,
                   "desired on-target interval size. Large on-target intervals will be split up into sub-interval of "
                   "about this size")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kAugmentBaitsDefaultOnTargetIntervalSize)
      ->group(group_name);
}

CLI::Option* AddOffTargetMinWidthOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding) {
  return app->add_option(opt::kOffTargetMinWidth, binding, "minimum width for off-targets")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kAugmentBaitsDefaultOffTargetMinWidth)
      ->group(group_name);
}

CLI::Option* AddOffTargetTrimLengthOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding) {
  return app
      ->add_option(opt::kOffTargetTrimLength,
                   binding,
                   "number of bases to trim off from candidate off-targets to prevent reads overlapping 2+ off-targets")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kAugmentBaitsDefaultOffTargetTrimLength)
      ->group(group_name);
}

CLI::Option* AddOffTargetIntervalSizeOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding) {
  return app
      ->add_option(opt::kOffTargetIntervalSize,
                   binding,
                   "desired off-target interval size. Large off-target intervals will be split up into sub-intervals "
                   "of about this size, while off-targets smaller than this size will be left as-is (if they exceed "
                   "off-target-min-width)")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kAugmentBaitsDefaultOffTargetIntervalSize)
      ->group(group_name);
}

CLI::Option* AddOnTargetMinMappabilityOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app->add_option(opt::kOnTargetMinMappability, binding, "minimum mappability threshold for on-target intervals")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(kAugmentBaitsDefaultOnTargetMinMappability)
      ->group(group_name);
}

CLI::Option* AddOffTargetMinMappabilityOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app
      ->add_option(opt::kOffTargetMinMappability, binding, "minimum mappability threshold for off-target intervals")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(kAugmentBaitsDefaultOffTargetMinMappability)
      ->group(group_name);
}

CLI::Option* AddSexChromosomeTargetMinMappabilityOption(const cli::AppPtr& app,
                                                        const std::string& group_name,
                                                        f64& binding) {
  return app
      ->add_option(opt::kSexChromosomeTargetMinMappability,
                   binding,
                   "minimum mappability threshold for sex-chromosome-target intervals. This is will be applied to "
                   "both on- and off-target intervals on sex chromosomes. Thus it overrides the "
                   "--on-target-min-mappability and --off-target-min-mappability values for targets on the sex "
                   "chromosome")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(kAugmentBaitsDefaultSexChromosomeTargetMinMappability)
      ->group(group_name);
}

// CalculateCoverage

CLI::Option* AddCoverageExcludeFlagsOption(const cli::AppPtr& app,
                                           const std::string& group_name,
                                           std::string& binding) {
  return app
      ->add_option(opt::kExcludeFlags,
                   binding,
                   "Exclude alignments if any bits in their FLAG field match the specified integer. Refer to "
                   "https://broadinstitute.github.io/picard/explain-flags.html for flag generation.")
      ->default_val(kCalculateCoverageDefaultExcludeFlags)
      ->group(group_name);
}

CLI::Option* AddCoverageIgnoreDnOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app
      ->add_flag(opt::kSkipDelAndRefBases, binding, "Whether to ignore deletion (D) and ref skip (N) in base counting.")
      ->default_val(kCalculateCoverageDefaultIgnoreDN)
      ->group(group_name);
}

// GCCorrect

CLI::Option* AddGcCorrectFirstSpanOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app
      ->add_option(opt::kFirstPassGcSmoothingFraction,
                   binding,
                   "Proportion of data points used for smoothing during the first LOESS iteration of GC bias "
                   "correction.")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(kGCCorrectDefaultFirstSpan)
      ->group(group_name);
}

// Denoise

CLI::Option* AddDenoiseMinTargetLengthOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding) {
  return app
      ->add_option(opt::kDenoiseMinTargetLength, binding, "minimum target size in panel to keep when calculating logrs")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kDenoiseDefaultMinTargetLength)
      ->group(group_name);
}

CLI::Option* AddDenoiseMinPanelMedianCoverageOption(const cli::AppPtr& app,
                                                    const std::string& group_name,
                                                    f64& binding) {
  return app
      ->add_option(opt::kDenoiseMinPanelMedianCoverage,
                   binding,
                   "filter away targets with minimum panel-of-normals median less than this value")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(kDenoiseDefaultMinPanelMedianCov)
      ->group(group_name);
}

CLI::Option* AddDenoiseMinPanelMedianAndTumorCoverageOption(const cli::AppPtr& app,
                                                            const std::string& group_name,
                                                            f64& binding) {
  return app
      ->add_option(opt::kDenoiseMinPanelMedianAndTumorCoverage,
                   binding,
                   "filter away targets with minimum total panel median + tumor count less than this value")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(kDenoiseDefaultMinPanelMedianAndTumorCov)
      ->group(group_name);
}

CLI::Option* AddDenoiseMinOffTargetFilterFractionOption(const cli::AppPtr& app,
                                                        const std::string& group_name,
                                                        f64& binding) {
  return app
      ->add_option(opt::kDenoiseMinOffTargetFilterFrac,
                   binding,
                   "if after filtering, the fraction of off-targets vs the total number of targets exceeds"
                   "this ratio then all off-targets will be ignored in logr output")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(kDenoiseDefaultMinOffTargetFilterFrac)
      ->group(group_name);
}

CLI::Option* AddDenoiseDisableFilterOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app->add_option(opt::kDenoiseDisableFilter, binding, "do not peform any filtering")
      ->check(CLI::TypeValidator<bool>())
      ->default_val(kDenoiseDefaultNoFilter)
      ->group(group_name);
}

CLI::Option* AddTumorCoverageOption(const cli::AppPtr& app,
                                    const std::string& group_name,
                                    std::optional<fs::path>& binding) {
  return app->add_option(opt::kTumorCoverage, binding, "tumor coverage")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddPanelOfNormalsCoveragesListOption(const cli::AppPtr& app,
                                                  const std::string& group_name,
                                                  fs::path& binding) {
  return app
      ->add_option(opt::kPanelOfNormalsCoveragesList, binding, "reference panel coverage file locations in a text file")
      ->check(cli::FileListReadableValidator())
      ->group(group_name);
}

// TwoSampleLogR / Shared coverage

CLI::Option* AddTumorMinCoverageOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding) {
  return app
      ->add_option(opt::kMinTumorCoverage,
                   binding,
                   "Minimum allowable coverage in the tumor sample for an interval to be kept. Intervals below "
                   "this value are excluded.")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kCopyNumberCallerDefaultTumorMinCoverage)
      ->group(group_name);
}

CLI::Option* AddNormalMinCoverageOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding) {
  return app
      ->add_option(opt::kMinNormalCoverage,
                   binding,
                   "Minimum allowable coverage in the normal sample for an interval to be kept. Intervals below "
                   "this value are excluded.")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kCopyNumberCallerDefaultSomaticNormalMinCoverage)
      ->group(group_name);
}

CLI::Option* AddTumorCoverageFileOption(const cli::AppPtr& app,
                                        const std::string& group_name,
                                        std::optional<fs::path>& binding) {
  return app->add_option(opt::kTumorCoverageFile, binding, "tumor coverage file")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddNormalCoverageFileOption(const cli::AppPtr& app,
                                         const std::string& group_name,
                                         std::optional<fs::path>& binding) {
  return app->add_option(opt::kNormalCoverageFile, binding, "normal coverage file")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddCoverageFileOption(const cli::AppPtr& app,
                                   const std::string& group_name,
                                   std::optional<fs::path>& binding) {
  return app->add_option(opt::kCoverageFile, binding, "coverage file")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

// Segmentation

CLI::Option* AddSeedSegmentsOption(const cli::AppPtr& app, const std::string& group_name, fs::path& binding) {
  return app
      ->add_option(opt::kSeedSegments,
                   binding,
                   "Path to pre-generated segments from which to start segmentation. The seed segments file must be in "
                   "SEG format.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddDisableHierarchicalPruningOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app
      ->add_flag(opt::kSkipHierarchicalPruning, binding, "Disable the post-segmentation hierarchical pruning step.")
      ->group(group_name);
}

CLI::Option* AddDisableMergingOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app
      ->add_flag(opt::kSkipTTestMerging,
                 binding,
                 "If enabled, the t-test comparison that merges sub-segments back together will be skipped. "
                 "Segments will not have a chance to be re-merged before continuing the recursive segmentation "
                 "process.")
      ->group(group_name);
}

CLI::Option* AddCbsMethodOption(const cli::AppPtr& app, const std::string& group_name, CbsMaxTMethod& binding) {
  return app
      ->add_option(opt::kCbsMethod,
                   binding,
                   "disable merging at the end of each recursive step. Only applicable if "
                   "--segmentation-mode=germline")
      ->transform(CLI::Transformer(segmentation::kStringToCbsMethod, CLI::ignore_case))
      ->group(group_name);
}

CLI::Option* AddSegmentationModeOption(const cli::AppPtr& app,
                                       const std::string& group_name,
                                       segmentation::SegmentationMode& binding) {
  return app->add_option(opt::kSegmentationMode, binding, "choose between germline or somatic")
      ->transform(CLI::Transformer(segmentation::kStringToSegmentationMode, CLI::ignore_case))
      ->group(group_name);
}

CLI::Option* AddTruncateOutliersOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app->add_flag(opt::kTruncateOutliers, binding, "disable merging at the end of each recursive step")
      ->group(group_name);
}

CLI::Option* AddPruningClusteringParameterOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app->add_option(opt::kHierarchicalPruningHeight, binding, "Clustering parameter for hierarchical pruning.")
      ->default_val(kSegmentationDefaultPruningClusteringParameter)
      ->group(group_name);
}

CLI::Option* AddMinTForAutomaticSegmentationOption(const cli::AppPtr& app,
                                                   const std::string& group_name,
                                                   f64& binding) {
  return app
      ->add_option(opt::kMinTValueToSkipPermutationTesting,
                   binding,
                   "Minimum t-test value threshold to skip significance testing for segmentation candidates.")
      ->default_val(kSegmentationDefaultMinTForAutomaticSegmentation)
      ->group(group_name);
}

CLI::Option* AddUndoLogrSegmentsOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app
      ->add_flag(opt::kEnableSegmentUndoing,
                 binding,
                 "Undo segments based on log-ratio values (i.e. re-merge segments with similar median log-ratio "
                 "values). This is done in an iterative fashion until `--max-segments-for-undoing` is reached or "
                 "no more segments can be undone based on the criteria set by `--initial-segment-undoing-factor`.")
      ->default_val(kSegmentationDefaultUndoLogrSegments)
      ->group(group_name);
}

CLI::Option* AddUndoLogrSegmentsSdFactorOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app
      ->add_option(opt::kInitialSegmentUndoingFactor,
                   binding,
                   "The initial standard deviation multiplier for segment undoing. Adjacent segments with a "
                   "difference less than or equal to the threshold determined by this multiplier are merged. "
                   "Only used when segment undoing is enabled.")
      ->default_val(kSegmentationDefaultUndoLogrSegmentsSdFactor)
      ->check(CLI::TypeValidator<f64>())
      ->check(CLI::PositiveNumber)
      ->group(group_name);
}

CLI::Option* AddIncrementUndoLogrSegmentsSdFactorOption(const cli::AppPtr& app,
                                                        const std::string& group_name,
                                                        f64& binding) {
  return app
      ->add_option(
          opt::kSegmentUndoingIncrement,
          binding,
          "The value added to `--initial-segment-undoing-factor` at each iteration. This gradually relaxes the "
          "merging criteria until the stop condition is met. Only used when segment undoing is enabled.")
      ->default_val(kSegmentationDefaultIncrementUndoLogrSegmentsSdFactor)
      ->check(CLI::TypeValidator<f64>())
      ->check(CLI::PositiveNumber)
      ->group(group_name);
}

CLI::Option* AddMaxNumSegmentsOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding) {
  return app
      ->add_option(opt::kMaxSegmentsForUndoing,
                   binding,
                   "Target segment count for the segment undoing procedure. The process stops once the total "
                   "number of segments is less than or equal to this value. "
                   "Only used when segment undoing is enabled.")
      ->default_val(kSegmentationDefaultMaxNumSegments)
      ->check(CLI::TypeValidator<size_t>())
      ->check(CLI::PositiveNumber)
      ->group(group_name);
}

CLI::Option* AddUndoDhSegmentsOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app
      ->add_flag(opt::kEnableDhSegmentUndoing,
                 binding,
                 "Undo segments based on DH values (re-merge segments with similar median DH values).")
      ->default_val(kSegmentationDefaultUndoDhSegments)
      ->group(group_name);
}

CLI::Option* AddUndoDhSegmentsSdFactorOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app
      ->add_option(opt::kInitialDhSegmentUndoingFactor,
                   binding,
                   "The initial standard deviation multiplier for DH segment undoing. Adjacent segments with a "
                   "difference less than or equal to the threshold determined by this multiplier are merged. "
                   "Requires `--enable-dh-segment-undoing`.")
      ->default_val(kSegmentationDefaultUndoDhSegmentsSdFactor)
      ->group(group_name);
}

CLI::Option* AddMinObservationsPerSegmentOption(const cli::AppPtr& app,
                                                const std::string& group_name,
                                                size_t& binding) {
  return app
      ->add_option(
          opt::kMinObservationsPerSegment, binding, "Minimum number of observations a segment is allowed to have")
      ->default_val(kSegmentationDefaultMinObsPerSegment)
      ->group(group_name);
}

// PurityPloidySearch

CLI::Option* AddSegMaxLogrOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app
      ->add_option(opt::kMaxSegmentLogRatioForPurityPloidySearch,
                   binding,
                   "Maximum log-ratio value allowed for a segment to be included in the purity/ploidy search. "
                   "Segments with a value greater than this are excluded.")
      ->default_val(kSegMaxLogR)
      ->group(group_name);
}

CLI::Option* AddSegMinNumLogrsOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding) {
  return app
      ->add_option(opt::kMinObservationsPerSegment,
                   binding,
                   "Minimum number of log-ratio observations required to include a segment in the analysis. "
                   "Segments with fewer observations than this value are excluded.")
      ->default_val(kSegMinNumLogRs)
      ->group(group_name);
}

CLI::Option* AddSegMinNumSnpsOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding) {
  return app
      ->add_option(opt::kMinSnpsPerSegmentForPurityPloidySearch,
                   binding,
                   "Minimum number of SNPs required to consider a segment in the purity/ploidy search. Segments "
                   "with fewer SNPs than this value are excluded.")
      ->default_val(kSegMinNumSnps)
      ->group(group_name);
}

// Sample Metadata

CLI::Option* AddSexOption(const cli::AppPtr& app, const std::string& group_name, std::optional<Sex>& binding) {
  return app->add_option(opt::kSex, binding, "Sex of sample.")
      ->transform(CLI::Transformer(kCharToSex, CLI::ignore_case))
      ->group(group_name);
}

CLI::Option* AddSampleIdOption(const cli::AppPtr& app, const std::string& group_name, std::string& binding) {
  return app->add_option(opt::kSampleId, binding, "Sample ID.")
      ->check(CLI::TypeValidator<std::string>())
      ->group(group_name);
}

CLI::Option* AddSampleNameOption(const cli::AppPtr& app,
                                 const std::string& group_name,
                                 std::optional<std::string>& binding) {
  return app
      ->add_option(
          opt::kSampleName,
          binding,
          "Name of sample in the input VCF header and used as a sample ID in the output. Does not require `--vcf`.")
      ->check(CLI::TypeValidator<std::string>())
      ->group(group_name);
}

CLI::Option* AddNormalSampleNameOption(const cli::AppPtr& app,
                                       const std::string& group_name,
                                       std::optional<std::string>& binding) {
  return app->add_option(opt::kNormalSampleName, binding, "Normal sample name.")
      ->check(CLI::TypeValidator<std::string>())
      ->group(group_name);
}

CLI::Option* AddTumorSampleNameOption(const cli::AppPtr& app,
                                      const std::string& group_name,
                                      std::optional<std::string>& binding) {
  return app->add_option(opt::kTumorSampleName, binding, "Tumor sample name.")
      ->check(CLI::TypeValidator<std::string>())
      ->group(group_name);
}

CLI::Option* AddTumorPurityOption(const cli::AppPtr& app, const std::string& group_name, std::optional<f64>& binding) {
  return app->add_option(opt::kTumorPurity, binding, "tumor purity")
      ->check(CLI::TypeValidator<f64>())
      ->check(CLI::Range(static_cast<f64>(0.15), static_cast<f64>(1)))
      ->group(group_name);
}

CLI::Option* AddTumorPloidyOption(const cli::AppPtr& app, const std::string& group_name, std::optional<f64>& binding) {
  return app->add_option(opt::kTumorPloidy, binding, "tumor ploidy")
      ->check(CLI::TypeValidator<f64>())
      ->check(CLI::PositiveNumber)
      ->group(group_name);
}

// Likelihood

CLI::Option* AddCnvLengthFlagMinSizeOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding) {
  return app->add_option(opt::kMinCnvLengthToFlag, binding, "Flag calls with a size lower than this threshold (2000).")
      ->default_val(kLikelihoodCnvLengthFlagMinSize)
      ->group(group_name);
}

// Save-flag helpers (new hidden advanced options)

CLI::Option* AddSaveTaskflowGraphOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app
      ->add_flag(opt::kSaveTaskflowGraph, binding, "Save intermediate taskflow graph output helpful for debugging.")
      ->group(group_name);
}

CLI::Option* AddSaveLogRatioSegmentsOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app
      ->add_flag(opt::kSaveLogRatioSegments,
                 binding,
                 "Output file for pre-likelihood LogR segments. These segments have gone through the undoing and "
                 "pruning processes already, but have not merged based on equal copy number or by mapping "
                 "qualities.")
      ->group(group_name);
}

CLI::Option* AddSaveBafSegmentsOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app
      ->add_flag(opt::kSaveBafSegments,
                 binding,
                 "Output file for segments after BAF segmentation but before hierarchical pruning.")
      ->group(group_name);
}

CLI::Option* AddMapqCutoffToFlagCallsOption(const cli::AppPtr& app, const std::string& group_name, s32& binding) {
  return app
      ->add_option(opt::kMinMapqForCalls,
                   binding,
                   "Minimum average mapping quality (MAPQ) required for a call in the SEG output. Segments with "
                   "an average MAPQ less than this value are merged. If merging is not possible, the call will be "
                   "flagged in the SEG output. Use a negative value to disable this filter.")
      ->default_val(kLikelihoodDefaultGermlineMAPQCutoffForCalls)
      ->group(group_name);
}

// MergeSegments

CLI::Option* AddMinMapqOption(const cli::AppPtr& app, const std::string& group_name, f64& binding) {
  return app
      ->add_option(
          opt::kMinMapq, binding, "minimum avg mean MAPQ for keeping segments. values less than this will get merged")
      ->default_val(kMergeDefaultMinMapq)
      ->group(group_name);
}

CLI::Option* AddUseMapqsFromFileOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app
      ->add_flag(opt::kUseMapqsFromFile,
                 binding,
                 "use the mean mapqs per target from the `--mapqs` file instead of the avg mean MAPQs per segment "
                 "contained the segments file. Can use this if the segment file doesn't contain a meanAvgMapq field")
      ->needs(opt::kMapqs)
      ->default_val(kMergeSegmentsDefaultUseMapqsObservations)
      ->group(group_name);
}

CLI::Option* AddRecalculatePerSegmentDataOption(const cli::AppPtr& app, const std::string& group_name, bool& binding) {
  return app
      ->add_flag(opt::kRecalculatePerSegmentData,
                 binding,
                 "recalculate the following fields from the SEG file: NumLogRatio, MeanLogRatio, NumSnp, MeanDh, "
                 "AvgMeanMapq etc. Otherwise, calculates the weighted mean of these fields for each merged segment.")
      ->needs(opt::kLogRatios)
      ->needs(opt::kMapqs)
      ->needs(opt::kSex)
      ->default_val(kMergeSegmentsDefaultRecalculatePerSegmentData)
      ->group(group_name);
}

CLI::Option* AddOutputDirOption(const cli::AppPtr& app, const std::string& group_name, fs::path& binding) {
  return app->add_option(opt::kOutputDir, binding, "Path to output directory.")->default_val(".")->group(group_name);
}

void DefineSomaticCommonOptions(const cli::AppPtr& app,
                                CopyNumberCallerOptions& options,
                                const SegmentUndoMode undo_mode) {
  const bool hide = true;
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options.output_dir);
  DefineIntervalOptions(app,
                        options.augment_baits_options,
                        hide,
                        kAugmentBaitsDefaultWholeGenomeIntervalSize,
                        kAugmentBaitsDefaultWholeGenomeMinMappability);
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options.sample_metadata_options.sex);
  AddTumorSampleNameOption(app, opt::kSampleMetadataOptionsGroup, options.sample_metadata_options.tumor_sample_name);
  cli::AddBamFileOption(app,
                        opt::kTumorBam,
                        options.tumor_bam_fname,
                        discarded_bam_index,
                        "Path to the input tumor BAM file. The BAM file must be coordinate-sorted and indexed.")
      ->required()
      ->check(cli::FileExtensionValidator({".bam"}))
      ->group(opt::kInputOptionsGroup);
  DefineAugmentBaitsInputOptions(app, options, !hide);
  AddSeedSegmentsOption(app, opt::kInputOptionsGroup, options.seed_segments_fname)->required();
  DefineCalculateCoverageAlgorithmOptions(app, options.calculate_coverage_options, hide);
  DefineGCCorrectAlgorithmOptions(app, options.gc_correct_options, hide);
  DefineSegmentationAlgorithmOptions(app, options.segmentation_options, hide, undo_mode);
  DefineSomaticSegmentationOptions(app, options.segmentation_options, hide);
  DefinePurityPloidySearchAlgorithmOptions(app, options.purity_ploidy_search_options, hide);
  DefineCoverageOptions(app, options);
  DefineBafFilteringOptions(app, options, hide);
  AddSaveTaskflowGraphOption(app, opt::kHiddenOptionsGroup, options.save_taskflow_graph);
  cli::AddThreadCountOption(app, opt::kThreads, options.threads)->group(opt::kPerformanceOptionsGroup);
}

void DefineSomaticTumorTargetedEnrichmentOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  DefineSomaticCommonOptions(app, *options, SegmentUndoMode::kOptIn);
  DefineDenoiseSubmoduleOptions(app, options->denoise_options);
  AddPanelOfNormalsOption(app, "SomaticTumorTargetedEnrichment", options->panel_of_normals_hdf5_fname)->required();
  AddIntervalsOption(app, "SomaticTumorTargetedEnrichment", options->augmented_baits_fname)->required();
  DefineBafFilteringOptions(app, *options, true);
}

void DefineSomaticTumorNormalWGSOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  cli::AddBamFileOption(app,
                        opt::kNormalBam,
                        options->normal_bam_fname,
                        discarded_bam_index,
                        "Path to the input normal BAM file. The BAM file must be coordinate-sorted and indexed.")
      ->required()
      ->check(cli::FileExtensionValidator({".bam"}))
      ->group(opt::kInputOptionsGroup);
  // Register --vcf before DefineSomaticCommonOptions so that VCF-dependent options (e.g. baf filter
  // options that ->needs("--vcf")) can find it. The required() and needs() constraints are applied
  // after --normal-sample-name and --tumor-sample-name are registered.
  auto* vcf_opt = AddVcfOption(app,
                               opt::kInputOptionsGroup,
                               options->vcf_fname,
                               "Path to an input VCF file with germline variant sites. Heterozygous germline calls "
                               "allow for estimation of allele-specific copy number.")
                      ->required();
  DefineSomaticCommonOptions(app, *options, SegmentUndoMode::kAlwaysEnabled);
  // Register --normal-sample-name after DefineSomaticCommonOptions (which registers --tumor-sample-name)
  // so vcf_opt can declare needs on both sample name options.
  AddNormalSampleNameOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.normal_sample_name);
  vcf_opt->needs(opt::kNormalSampleName)->needs(opt::kTumorSampleName);
  // We also include augment-baits options here because tumor-normal will require an
  // augmented-baits file created from the whole genome. We do not require the user to run AugmentBaits themselves for
  // tumor-normal. (for TE, AugmentBaits is run as a part of GeneratePanelOfNormals)
  AddIntervalsOption(app, opt::kHiddenOptionsGroup, options->augmented_baits_fname);
  AddSaveBafSegmentsOption(app, opt::kHiddenOptionsGroup, options->save_baf_segments);
  // Register --tumor-ploidy before --tumor-purity, then wire up mutual needs after both are registered
  auto* ploidy_opt =
      AddTumorPloidyOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.ploidy);
  auto* purity_opt =
      AddTumorPurityOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.purity);
  purity_opt->needs(ploidy_opt);
  ploidy_opt->needs(purity_opt);
}

void DefineGermlineNormalWGSOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  // hide options for advanced users
  const bool hide = true;
  cli::AddBamFileOption(app,
                        opt::kBam,
                        options->normal_bam_fname,
                        discarded_bam_index,
                        "Path to the input BAM file. The BAM file must be coordinate-sorted and indexed.")
      ->required()
      ->check(cli::FileExtensionValidator({".bam"}))
      ->group(opt::kInputOptionsGroup);
  DefineAugmentBaitsInputOptions(app, *options, false);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  DefineIntervalOptions(app,
                        options->augment_baits_options,
                        hide,
                        kAugmentBaitsGermlineWGSDefaultIntervalSize,
                        kAugmentBaitsGermlineWGSDefaultWholeGenomeMinMappabiilty);
  AddSeedSegmentsOption(app, opt::kInputOptionsGroup, options->seed_segments_fname)->required();
  // Register --sample-name before --vcf so VCF can declare a dependency on it
  AddSampleNameOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.normal_sample_name);
  AddVcfOption(app,
               opt::kInputOptionsGroup,
               options->vcf_fname,
               "Path to an input VCF file for germline SNPs (for BAF visualization).")
      ->needs(opt::kSampleName);
  DefineCalculateCoverageAlgorithmOptions(app, options->calculate_coverage_options, hide);
  DefineSegmentationAlgorithmOptions(app, options->segmentation_options, hide, SegmentUndoMode::kOptIn);
  DefineGermlineCallFlagOptions(app, options->likelihood_options, hide);
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex);
  cli::AddThreadCountOption(app, opt::kThreads, options->threads)->group(opt::kPerformanceOptionsGroup);
  // Hidden advanced options
  AddNormalSampleMinDepthOption(app,
                                opt::kHiddenOptionsGroup,
                                options->baf_filter_options.normal_sample_min_depth,
                                kVcfParsingOptionsDefaultGermlineNormalSampleMinDepth);
  AddGcCorrectFirstSpanOption(app, opt::kHiddenOptionsGroup, options->gc_correct_options.first_span);
  AddIntervalsOption(app, opt::kHiddenOptionsGroup, options->augmented_baits_fname);
  AddSaveLogRatioSegmentsOption(app, opt::kHiddenOptionsGroup, options->save_log_ratio_segments);
  AddSaveTaskflowGraphOption(app, opt::kHiddenOptionsGroup, options->save_taskflow_graph);
}

// ==============================================================================
// AugmentBaits
// ==============================================================================

void DefineAugmentBaitsInputOptions(const cli::AppPtr& app, CopyNumberCallerOptions& options, const bool hide) {
  cli::AddFastaFileOption(app,
                          opt::kReference,
                          options.reference_genome_fname,
                          "Path to the input reference FASTA file. The FASTA file must be indexed with a .fai file.")
      ->group(opt::kInputOptionsGroup)
      ->required();
  AddMappabilityBigwigOption(app, opt::kInputOptionsGroup, options.mappability_bigwig_fname)->required();
  AddBlocklistBedOption(app, hide ? opt::kHiddenOptionsGroup : opt::kInputOptionsGroup, options.blocklist_bed_fname);
}

void DefineIntervalOptions(const cli::AppPtr& app,
                           AugmentBaitsOptions& options,
                           const bool hide,
                           const size_t default_interval_size,
                           const f64 default_mappability) {
  AddWholeGenomeIntervalSizeOption(
      app, opt::kRegionOptionsGroup, options.whole_genome_interval_size, default_interval_size);
  AddWholeGenomeMinMappabilityOption(
      app, opt::kRegionOptionsGroup, options.whole_genome_min_mappability, default_mappability);
  AddMinGcContentOption(app, hide ? opt::kHiddenOptionsGroup : opt::kRegionOptionsGroup, options.min_gc_content);
  AddMaxGcContentOption(app, hide ? opt::kHiddenOptionsGroup : opt::kRegionOptionsGroup, options.max_gc_content);
}

void DefineAdvancedIntervalOptions(const cli::AppPtr& app, AugmentBaitsOptions& options) {
  AddNoOffTargetsOption(app, opt::kRegionOptionsGroup, options.no_off_targets);
  AddOnTargetIntervalSizeOption(app, opt::kRegionOptionsGroup, options.on_target_interval_size);
  AddOffTargetMinWidthOption(app, opt::kRegionOptionsGroup, options.off_target_min_width);
  AddOffTargetTrimLengthOption(app, opt::kRegionOptionsGroup, options.off_target_trim_length);
  AddOffTargetIntervalSizeOption(app, opt::kRegionOptionsGroup, options.off_target_interval_size);
  AddOnTargetMinMappabilityOption(app, opt::kRegionOptionsGroup, options.on_target_min_mappability);
  AddOffTargetMinMappabilityOption(app, opt::kRegionOptionsGroup, options.off_target_min_mappability);
  AddSexChromosomeTargetMinMappabilityOption(app, opt::kRegionOptionsGroup, options.sex_chromosome_min_mappability);
}

void DefineAugmentBaitsSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  DefineAugmentBaitsInputOptions(app, *options, false);
  cli::AddBamFileOption(app,
                        opt::kBam,
                        options->normal_bam_fname,
                        discarded_bam_index,
                        "Path to the input BAM file. Used to check whether the chrY PAR region was masked during "
                        "alignment. If unmasked, PAR intervals are removed from the output baits file.")
      ->required()
      ->check(cli::FileExtensionValidator({".bam"}))
      ->group(opt::kInputOptionsGroup);
  AddSeedSegmentsOption(app, opt::kInputOptionsGroup, options->seed_segments_fname)->required();
  AddBedOption(app, opt::kInputOptionsGroup, options->augmented_baits_fname);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddWholeGenomeOption(app, opt::kRegionOptionsGroup, options->augment_baits_options.generate_whole_genome_baits);
  DefineIntervalOptions(app,
                        options->augment_baits_options,
                        false,
                        kAugmentBaitsDefaultWholeGenomeIntervalSize,
                        kAugmentBaitsDefaultWholeGenomeMinMappability);
  DefineAdvancedIntervalOptions(app, options->augment_baits_options);
  AddIncludeAltContigsOption(app, opt::kRegionOptionsGroup, options->augment_baits_options.include_alt_contigs);
}

// ==============================================================================
// CalculateCoverage
// ==============================================================================

void DefineCalculateCoverageAlgorithmOptions(const cli::AppPtr& app,
                                             CalculateCoverageOptions& options,
                                             const bool hide) {
  // we are not validating this as an s32 yet; currently validation happens with an htslib call down stream
  AddCoverageExcludeFlagsOption(app, opt::kReadFilteringOptionsGroup, options.exclude_flags);
  AddCoverageIgnoreDnOption(
      app, hide ? opt::kHiddenOptionsGroup : opt::kCalculateCoverageOptionsGroup, options.ignore_DN);
}

void DefineCalculateCoverageSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  cli::AddBamFileOption(app,
                        opt::kNormalBam,
                        options->normal_bam_fname,
                        discarded_bam_index,
                        "Path to the input BAM file. The BAM file must be coordinate-sorted and indexed.")
      ->required()
      ->check(cli::FileExtensionValidator({".bam"}))
      ->group(opt::kInputOptionsGroup);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddIntervalsOption(app, opt::kInputOptionsGroup, options->augmented_baits_fname)->required();
  DefineCalculateCoverageAlgorithmOptions(app, options->calculate_coverage_options, false);
  cli::AddThreadCountOption(app, opt::kThreads, options->threads)->group(opt::kPerformanceOptionsGroup);
}

// ==============================================================================
// GCCorrect
// ==============================================================================

void DefineGCCorrectAlgorithmOptions(const cli::AppPtr& app, GCCorrectOptions& options, const bool hide) {
  const auto group_name = hide ? opt::kHiddenOptionsGroup : opt::kGCCorrectOptionsGroup;
  AddGcCorrectFirstSpanOption(app, group_name, options.first_span);
}

void DefineGCCorrectSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  AddCoverageFileOption(app, opt::kInputOptionsGroup, options->normal_coverage_fname)->required();
  AddIntervalsOption(app, opt::kInputOptionsGroup, options->augmented_baits_fname)->required();
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  DefineGCCorrectAlgorithmOptions(app, options->gc_correct_options, false);
  cli::AddThreadCountOption(app, opt::kThreads, options->threads)->group(opt::kPerformanceOptionsGroup);
}

// ==============================================================================
// Denoise
// ==============================================================================

void DefineDenoiseSubmoduleOptions(const cli::AppPtr& app, DenoiseOptions& options) {
  AddDenoiseMinTargetLengthOption(app, opt::kDenoiseOptionsGroup, options.min_target_length);
  AddDenoiseMinPanelMedianCoverageOption(app, opt::kDenoiseOptionsGroup, options.min_panel_median_cov);
  AddDenoiseMinPanelMedianAndTumorCoverageOption(
      app, opt::kDenoiseOptionsGroup, options.min_panel_median_and_tumor_cov);
  AddDenoiseMinOffTargetFilterFractionOption(app, opt::kDenoiseOptionsGroup, options.min_off_target_filter_frac);
  AddDenoiseDisableFilterOption(app, opt::kDenoiseOptionsGroup, options.no_filter);
}

void DefineDenoiseSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  AddTumorCoverageOption(app, opt::kInputOptionsGroup, options->tumor_coverage_fname)->required();
  AddPanelOfNormalsCoveragesListOption(app, opt::kInputOptionsGroup, options->panel_of_normals_lists)->required();
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  DefineDenoiseSubmoduleOptions(app, options->denoise_options);
}

// ==============================================================================
// TwoSampleLogR
// ==============================================================================

void DefineCoverageOptions(const cli::AppPtr& app, CopyNumberCallerOptions& options) {
  AddTumorMinCoverageOption(app, opt::kReadFilteringOptionsGroup, options.tumor_min_coverage);
  AddNormalMinCoverageOption(app, opt::kReadFilteringOptionsGroup, options.normal_min_coverage);
}

void DefineTwoSampleLogRSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  AddTumorCoverageFileOption(app, opt::kInputOptionsGroup, options->tumor_coverage_fname)->required();
  AddNormalCoverageFileOption(app, opt::kInputOptionsGroup, options->normal_coverage_fname)->required();
  cli::AddFastaFileOption(app,
                          opt::kReference,
                          options->reference_genome_fname,
                          "Path to the input reference FASTA file. The FASTA file must be indexed with a .fai file.")
      ->group(opt::kInputOptionsGroup);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  DefineCoverageOptions(app, *options);
}

// ==============================================================================
// OneSampleLogR
// ==============================================================================

void DefineOneSampleLogRSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  AddCoverageFileOption(app, opt::kInputOptionsGroup, options->normal_coverage_fname)->required();
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
}

// ==============================================================================
// Segmentation
// ==============================================================================

void DefineSegmentationAlgorithmOptions(const cli::AppPtr& app,
                                        SegmentationOptions& options,
                                        const bool hide,
                                        const SegmentUndoMode undo_mode) {
  const auto group_name = hide ? opt::kHiddenOptionsGroup : opt::kSegmentationOptionsGroup;
  AddPruningClusteringParameterOption(app, opt::kSegmentationOptionsGroup, options.pruning_clustering_parameter);
  if (undo_mode == SegmentUndoMode::kOptIn) {
    AddUndoLogrSegmentsOption(app, opt::kSegmentationOptionsGroup, options.undo_logr_segments);
  }
  AddUndoLogrSegmentsSdFactorOption(app, opt::kSegmentationOptionsGroup, options.undo_logr_segments_sd_factor);
  AddMaxNumSegmentsOption(app, opt::kSegmentationOptionsGroup, options.max_num_segments);
  AddIncrementUndoLogrSegmentsSdFactorOption(app, group_name, options.increment_undo_logr_segments_sd_factor);
  AddDisableMergingOption(app, group_name, options.disable_merging);
  AddDisableHierarchicalPruningOption(app, group_name, options.disable_hierarchical_pruning);
}

void DefineSomaticSegmentationOptions(const cli::AppPtr& app, SegmentationOptions& options, const bool hide) {
  const auto group_name = hide ? opt::kHiddenOptionsGroup : opt::kSegmentationOptionsGroup;
  AddMinTForAutomaticSegmentationOption(app, group_name, options.min_t_for_automatic_segmentation);
  AddUndoDhSegmentsOption(app, group_name, options.undo_dh_segments);
  AddUndoDhSegmentsSdFactorOption(app, group_name, options.undo_dh_segments_sd_factor)
      ->needs(opt::kEnableDhSegmentUndoing);
}

void DefineBafFilteringOptions(const cli::AppPtr& app, CopyNumberCallerOptions& options, const bool hide) {
  const auto group_name = hide ? opt::kHiddenOptionsGroup : opt::kVcfParsingOptionsGroup;
  AddNormalSampleMinDepthOption(app,
                                opt::kVcfParsingOptionsGroup,
                                options.baf_filter_options.normal_sample_min_depth,
                                kVcfParsingOptionsDefaultSomaticNormalSampleMinDepth);
  AddNormalSampleMinBafOption(app, opt::kVcfParsingOptionsGroup, options.baf_filter_options.normal_sample_min_baf);
  AddNormalSampleMaxBafOption(app, opt::kVcfParsingOptionsGroup, options.baf_filter_options.normal_sample_max_baf);
  AddTumorSampleMinDepthOption(app, group_name, options.baf_filter_options.tumor_sample_min_depth);
  AddTumorSampleMinBafOption(app, group_name, options.baf_filter_options.tumor_sample_min_baf);
  AddTumorSampleMaxBafOption(app, group_name, options.baf_filter_options.tumor_sample_max_baf);
  AddAdInOption(app, group_name, options.ad_fname);
  AddForceEnableSomaticVariantParsingOption(
      app, opt::kVcfParsingOptionsGroup, options.baf_filter_options.force_enable_somatic_variant_parsing);
}

void DefineSegmentationSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  AddLogrsOption(app, opt::kInputOptionsGroup, options->logrs_fname)->required();
  AddSeedSegmentsOption(app, opt::kInputOptionsGroup, options->seed_segments_fname)->required();
  // Register --vcf before --reference so the fasta option can declare ->needs("--vcf")
  AddVcfOption(
      app, opt::kInputOptionsGroup, options->vcf_fname, "Path to an input VCF file with germline variant sites.");
  cli::AddFastaFileOption(app,
                          opt::kReference,
                          options->reference_genome_fname,
                          "Path to the input reference FASTA file. The FASTA file must be indexed with a .fai file.")
      ->group(opt::kInputOptionsGroup)
      ->needs(opt::kVcf);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddSegmentationModeOption(app, opt::kSegmentationMode, options->segmentation_options.segmentation_mode)->required();
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex)->required();
  AddNormalSampleNameOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.normal_sample_name);
  AddTumorSampleNameOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.tumor_sample_name);
  AddMinObservationsPerSegmentOption(
      app, opt::kSegmentationOptionsGroup, options->segmentation_options.min_obs_per_segment);
  DefineSegmentationAlgorithmOptions(app, options->segmentation_options, false, SegmentUndoMode::kOptIn);
  DefineSomaticSegmentationOptions(app, options->segmentation_options, false);
  cli::AddThreadCountOption(app, opt::kThreads, options->threads)->group(opt::kPerformanceOptionsGroup);
  DefineBafFilteringOptions(app, *options, false);
  AddCbsMethodOption(app, opt::kSegmentationOptionsGroup, options->segmentation_options.cbs_method);
  AddTruncateOutliersOption(app, opt::kSegmentationOptionsGroup, options->segmentation_options.truncate_outliers);
}

// ==============================================================================
// PurityPloidySearch
// ==============================================================================

void DefinePurityPloidySearchAlgorithmOptions(const cli::AppPtr& app,
                                              PurityPloidySearchOptions& options,
                                              const bool hide) {
  const auto group_name = hide ? opt::kHiddenOptionsGroup : opt::kPurityPloidySearchOptionsGroup;
  AddSegMaxLogrOption(app, group_name, options.seg_max_logr);
  AddSegMinNumLogrsOption(app, group_name, options.seg_min_num_logrs);
  AddSegMinNumSnpsOption(app, group_name, options.seg_min_num_snps);
}

void DefinePurityPloidySearchSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  AddSegmentsOption(app, opt::kInputOptionsGroup, options->segments_fname)->required();
  AddLogrsOption(app, opt::kInputOptionsGroup, options->logrs_fname)->required();
  AddBafsOption(app, opt::kInputOptionsGroup, options->bafs_fname)->required();
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex)->required();
  AddSampleIdOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sample_id);
  DefinePurityPloidySearchAlgorithmOptions(app, options->purity_ploidy_search_options, false);
}

// ==============================================================================
// Likelihood (PredictSomaticCNA / PredictGermlineCNV)
// ==============================================================================

void DefinePredictSomaticCNASubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  AddSegmentsOption(app, opt::kInputOptionsGroup, options->segments_fname)->required();
  AddLogrsOption(app, opt::kInputOptionsGroup, options->logrs_fname)->required();
  AddBafsOption(app, opt::kInputOptionsGroup, options->bafs_fname)->required();
  cli::AddFastaFileOption(app,
                          opt::kReference,
                          options->reference_genome_fname,
                          "Path to the input reference FASTA file. The FASTA file must be indexed with a .fai file.")
      ->group(opt::kInputOptionsGroup)
      ->required();
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  // Register --tumor-ploidy before --tumor-purity, then wire up mutual needs after both are registered
  auto* ploidy_opt =
      AddTumorPloidyOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.ploidy)->required();
  auto* purity_opt =
      AddTumorPurityOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.purity)->required();
  purity_opt->needs(ploidy_opt);
  ploidy_opt->needs(purity_opt);
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex)->required();
}

void DefineGermlineCallFlagOptions(const cli::AppPtr& app, LikelihoodOptions& options, const bool hide) {
  AddMapqCutoffToFlagCallsOption(app, opt::kOutputOptionsGroup, options.mapq_cutoff_for_calls);
  AddCnvLengthFlagMinSizeOption(
      app, hide ? opt::kHiddenOptionsGroup : opt::kOutputOptionsGroup, options.cnv_length_flag_min_size);
}

void DefinePredictGermlineCNVSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  AddSegmentsOption(app, opt::kInputOptionsGroup, options->segments_fname)->required();
  AddLogrsOption(app, opt::kInputOptionsGroup, options->logrs_fname)->required();
  AddMapqsOption(app, opt::kInputOptionsGroup, options->mapping_qualities_fname);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex)->required();
  AddSampleIdOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sample_id);
  DefineGermlineCallFlagOptions(app, options->likelihood_options, false);
}

// ==============================================================================
// MergeSegments
// ==============================================================================

void DefineMergeSegmentsSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  AddSegmentsOption(app, opt::kInputOptionsGroup, options->segments_fname)->required();
  AddLogrsOption(app, opt::kInputOptionsGroup, options->logrs_fname)->required();
  AddMapqsOption(app, opt::kInputOptionsGroup, options->mapping_qualities_fname)->required();
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddUseMapqsFromFileOption(app, opt::kInputOptionsGroup, options->merge_segments_options.use_mapqs_observations);
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex)->required();
  AddSampleIdOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sample_id);
  AddMinMapqOption(app, opt::kOutputOptionsGroup, options->merge_segments_options.min_mapq_threshold);
  AddRecalculatePerSegmentDataOption(
      app, opt::kOutputOptionsGroup, options->merge_segments_options.recalculate_per_segment_data);
}

// ==============================================================================
// SegToVcf
// ==============================================================================

void DefineSegToVcfSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options) {
  AddSegmentsOption(app, opt::kInputOptionsGroup, options->segments_fname)->required();
  cli::AddFastaFileOption(app,
                          opt::kReference,
                          options->reference_genome_fname,
                          "Path to the input reference FASTA file. The FASTA file must be indexed with a .fai file.")
      ->group(opt::kInputOptionsGroup);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex)->required();
  AddSampleIdOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sample_id);
}

}  // namespace xoos::cnc
