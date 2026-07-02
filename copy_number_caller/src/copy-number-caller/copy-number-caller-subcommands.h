#pragma once

#include <xoos/cli/cli.h>
#include <xoos/types/int.h>

namespace xoos::cnc {

constexpr auto* kSomaticTumorNormalWGSSubcommand = "tumor-normal-wgs";
constexpr auto* kGermlineNormalWGSSubcommand = "germline-wgs";
constexpr auto* kAugmentBaitsSubcommand = "augment-baits";
constexpr auto* kCalculateCoverageSubcommand = "calculate-coverage";
constexpr auto* kGCCorrectSubcommand = "gc-correct";
constexpr auto* kPredictSomaticCNASubcommand = "predict-somatic-cnv";
constexpr auto* kPredictGermlineCNVSubcommand = "predict-germline-cnv";
constexpr auto* kDenoiseSubcommand = "denoise";
constexpr auto* kTwoSampleLogRSubcommand = "two-sample-logr";
constexpr auto* kOneSampleLogRSubcommand = "one-sample-logr";
constexpr auto* kSegmentationSubcommand = "segmentation";
constexpr auto* kPurityPloidySearchSubcommand = "purity-ploidy-search";
constexpr auto* kMergeSegmentsSubcommand = "merge-segments";
constexpr auto* kSegToVcfSubcommand = "seg-to-vcf";

}  // namespace xoos::cnc
