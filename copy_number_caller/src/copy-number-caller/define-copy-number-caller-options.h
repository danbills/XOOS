#pragma once

#include <string>

#include <CLI/App.hpp>

#include <xoos/cli/cli.h>

#include "copy-number-caller/copy-number-caller-options.h"  // IWYU pragma: keep

namespace xoos::cnc {

// End-to-end pipeline options
void DefineSomaticTumorTargetedEnrichmentOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);
void DefineSomaticTumorNormalWGSOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);
void DefineGermlineNormalWGSOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);

// PredictSomaticCNA options
void DefinePredictSomaticCNASubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);

// Controls whether --enable-segment-undoing is registered as a CLI flag.
// kAlwaysEnabled: segment undoing is unconditionally on; no flag is registered.
// kOptIn: the flag is registered so users can opt in (default: off).
enum class SegmentUndoMode {
  kOptIn,
  kAlwaysEnabled
};

// Shared options
void DefineSomaticCommonOptions(const cli::AppPtr& app, CopyNumberCallerOptions& options, SegmentUndoMode undo_mode);
void DefineSomaticSegmentationOptions(const cli::AppPtr& app, SegmentationOptions& options, bool hide);
void DefineSegmentationAlgorithmOptions(const cli::AppPtr& app,
                                        SegmentationOptions& options,
                                        bool hide,
                                        SegmentUndoMode undo_mode);
void DefineBafFilteringOptions(const cli::AppPtr& app, CopyNumberCallerOptions& options, bool hide);
void DefineGermlineCallFlagOptions(const cli::AppPtr& app, LikelihoodOptions& options, bool hide);
void DefineAugmentBaitsInputOptions(const cli::AppPtr& app, CopyNumberCallerOptions& options, bool hide);
void DefineIntervalOptions(const cli::AppPtr& app,
                           AugmentBaitsOptions& options,
                           bool hide,
                           size_t default_interval_size,
                           f64 default_mappability);
void DefineAdvancedIntervalOptions(const cli::AppPtr& app, AugmentBaitsOptions& options);
void DefineCalculateCoverageAlgorithmOptions(const cli::AppPtr& app, CalculateCoverageOptions& options, bool hide);
void DefineGCCorrectAlgorithmOptions(const cli::AppPtr& app, GCCorrectOptions& options, bool hide);
void DefineDenoiseSubmoduleOptions(const cli::AppPtr& app, DenoiseOptions& options);
void DefineCoverageOptions(const cli::AppPtr& app, CopyNumberCallerOptions& options);
void DefinePurityPloidySearchAlgorithmOptions(const cli::AppPtr& app, PurityPloidySearchOptions& options, bool hide);

// Hidden subcommand options
void DefineSegmentationSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);
void DefinePredictGermlineCNVSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);
void DefineAugmentBaitsSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);
void DefineCalculateCoverageSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);
void DefineGCCorrectSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);
void DefineDenoiseSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);
void DefineOneSampleLogRSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);
void DefineTwoSampleLogRSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);
void DefinePurityPloidySearchSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);
void DefineMergeSegmentsSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);
void DefineSegToVcfSubcommandOptions(const cli::AppPtr& app, CopyNumberCallerOptionsPtr& options);

// =============================================================================
// Add*Option helper functions
// =============================================================================

// Shared End-to-End options
CLI::Option* AddPanelOfNormalsOption(const cli::AppPtr& app, const std::string& group_name, fs::path& binding);
CLI::Option* AddIntervalsOption(const cli::AppPtr& app,
                                const std::string& group_name,
                                std::optional<fs::path>& binding);

// Shared IO options
CLI::Option* AddSegmentsOption(const cli::AppPtr& app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddLogrsOption(const cli::AppPtr& app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddBafsOption(const cli::AppPtr& app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddMapqsOption(const cli::AppPtr& app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddVcfOption(const cli::AppPtr& app,
                          const std::string& group_name,
                          std::optional<fs::path>& binding,
                          const std::string& description);
CLI::Option* AddOutputDirOption(const cli::AppPtr& app, const std::string& group_name, fs::path& binding);

// VCF Parsing
CLI::Option* AddNormalSampleMinDepthOption(const cli::AppPtr& app,
                                           const std::string& group_name,
                                           s32& binding,
                                           s32 default_depth);
CLI::Option* AddNormalSampleMinBafOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddNormalSampleMaxBafOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddTumorSampleMinDepthOption(const cli::AppPtr& app, const std::string& group_name, s32& binding);
CLI::Option* AddTumorSampleMinBafOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddTumorSampleMaxBafOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddAdInOption(const cli::AppPtr& app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddForceEnableSomaticVariantParsingOption(const cli::AppPtr& app,
                                                       const std::string& group_name,
                                                       bool& binding);

// AugmentBaits Inputs
CLI::Option* AddMappabilityBigwigOption(const cli::AppPtr& app, const std::string& group_name, fs::path& binding);
CLI::Option* AddBlocklistBedOption(const cli::AppPtr& app,
                                   const std::string& group_name,
                                   std::optional<fs::path>& binding);

// AugmentBaits Other
CLI::Option* AddIncludeAltContigsOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);

// AugmentBaits Whole Genome
CLI::Option* AddWholeGenomeOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);
CLI::Option* AddWholeGenomeIntervalSizeOption(const cli::AppPtr& app,
                                              const std::string& group_name,
                                              size_t& binding,
                                              size_t default_size);
CLI::Option* AddMinGcContentOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddMaxGcContentOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddWholeGenomeMinMappabilityOption(const cli::AppPtr& app,
                                                const std::string& group_name,
                                                f64& binding,
                                                f64 default_mappability);

// AugmentBaits Baits Options
CLI::Option* AddBedOption(const cli::AppPtr& app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddNoOffTargetsOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);
CLI::Option* AddOnTargetIntervalSizeOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding);
CLI::Option* AddOffTargetMinWidthOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding);
CLI::Option* AddOffTargetTrimLengthOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding);
CLI::Option* AddOffTargetIntervalSizeOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding);
CLI::Option* AddOnTargetMinMappabilityOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddOffTargetMinMappabilityOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddSexChromosomeTargetMinMappabilityOption(const cli::AppPtr& app,
                                                        const std::string& group_name,
                                                        f64& binding);

// CalculateCoverage
CLI::Option* AddCoverageExcludeFlagsOption(const cli::AppPtr& app, const std::string& group_name, std::string& binding);
CLI::Option* AddCoverageIgnoreDnOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);

// GCCorrect
CLI::Option* AddGcCorrectFirstSpanOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);

// Denoise
CLI::Option* AddDenoiseMinTargetLengthOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding);
CLI::Option* AddDenoiseMinPanelMedianCoverageOption(const cli::AppPtr& app,
                                                    const std::string& group_name,
                                                    f64& binding);
CLI::Option* AddDenoiseMinPanelMedianAndTumorCoverageOption(const cli::AppPtr& app,
                                                            const std::string& group_name,
                                                            f64& binding);
CLI::Option* AddDenoiseMinOffTargetFilterFractionOption(const cli::AppPtr& app,
                                                        const std::string& group_name,
                                                        f64& binding);
CLI::Option* AddDenoiseDisableFilterOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);
CLI::Option* AddTumorCoverageOption(const cli::AppPtr& app,
                                    const std::string& group_name,
                                    std::optional<fs::path>& binding);
CLI::Option* AddPanelOfNormalsCoveragesListOption(const cli::AppPtr& app,
                                                  const std::string& group_name,
                                                  fs::path& binding);

// TwoSampleLogR / Shared coverage
CLI::Option* AddTumorMinCoverageOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding);
CLI::Option* AddNormalMinCoverageOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding);
CLI::Option* AddTumorCoverageFileOption(const cli::AppPtr& app,
                                        const std::string& group_name,
                                        std::optional<fs::path>& binding);
CLI::Option* AddNormalCoverageFileOption(const cli::AppPtr& app,
                                         const std::string& group_name,
                                         std::optional<fs::path>& binding);
CLI::Option* AddCoverageFileOption(const cli::AppPtr& app,
                                   const std::string& group_name,
                                   std::optional<fs::path>& binding);

// Segmentation
CLI::Option* AddSeedSegmentsOption(const cli::AppPtr& app, const std::string& group_name, fs::path& binding);
CLI::Option* AddDisableHierarchicalPruningOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);
CLI::Option* AddDisableMergingOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);
CLI::Option* AddCbsMethodOption(const cli::AppPtr& app, const std::string& group_name, CbsMaxTMethod& binding);
CLI::Option* AddSegmentationModeOption(const cli::AppPtr& app,
                                       const std::string& group_name,
                                       segmentation::SegmentationMode& binding);
CLI::Option* AddTruncateOutliersOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);
CLI::Option* AddPruningClusteringParameterOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddMinTForAutomaticSegmentationOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddUndoLogrSegmentsOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);
CLI::Option* AddUndoLogrSegmentsSdFactorOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddIncrementUndoLogrSegmentsSdFactorOption(const cli::AppPtr& app,
                                                        const std::string& group_name,
                                                        f64& binding);
CLI::Option* AddMaxNumSegmentsOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding);
CLI::Option* AddUndoDhSegmentsOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);
CLI::Option* AddUndoDhSegmentsSdFactorOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddFfpeModeOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);
CLI::Option* AddMinObservationsPerSegmentOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding);

// PurityPloidySearch
CLI::Option* AddSegMaxLogrOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddSegMinNumLogrsOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding);
CLI::Option* AddSegMinNumSnpsOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding);

// Sample Metadata
CLI::Option* AddSexOption(const cli::AppPtr& app, const std::string& group_name, std::optional<Sex>& binding);
CLI::Option* AddSampleNameOption(const cli::AppPtr& app,
                                 const std::string& group_name,
                                 std::optional<std::string>& binding);
CLI::Option* AddTumorPurityOption(const cli::AppPtr& app, const std::string& group_name, std::optional<f64>& binding);
CLI::Option* AddTumorPloidyOption(const cli::AppPtr& app, const std::string& group_name, std::optional<f64>& binding);
CLI::Option* AddSampleIdOption(const cli::AppPtr& app, const std::string& group_name, std::string& binding);
CLI::Option* AddNormalSampleNameOption(const cli::AppPtr& app,
                                       const std::string& group_name,
                                       std::optional<std::string>& binding);
CLI::Option* AddTumorSampleNameOption(const cli::AppPtr& app,
                                      const std::string& group_name,
                                      std::optional<std::string>& binding);

// Likelihood
CLI::Option* AddCnvLengthFlagMinSizeOption(const cli::AppPtr& app, const std::string& group_name, size_t& binding);
CLI::Option* AddMapqCutoffToFlagCallsOption(const cli::AppPtr& app, const std::string& group_name, s32& binding);

// Save-flag helpers (new hidden advanced options)
CLI::Option* AddSaveTaskflowGraphOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);
CLI::Option* AddSaveLogRatioSegmentsOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);
CLI::Option* AddSaveBafSegmentsOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);

// MergeSegments
CLI::Option* AddMinMapqOption(const cli::AppPtr& app, const std::string& group_name, f64& binding);
CLI::Option* AddUseMapqsFromFileOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);
CLI::Option* AddRecalculatePerSegmentDataOption(const cli::AppPtr& app, const std::string& group_name, bool& binding);

}  // namespace xoos::cnc
