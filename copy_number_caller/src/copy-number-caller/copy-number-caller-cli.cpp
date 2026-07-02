#include "copy-number-caller/copy-number-caller-cli.h"

#include <algorithm>
#include <array>

#include <xoos/log/logging.h>
#include <xoos/util/file-functions.h>

#include "augment-baits/augment-baits.h"
#include "calculate-coverage/calculate-coverage.h"
#include "copy-number-caller/copy-number-caller-cli-option-names.h"
#include "copy-number-caller/copy-number-caller-subcommands.h"
#include "copy-number-caller/define-copy-number-caller-options.h"
#include "denoise/denoise.h"
#include "gc-correct/gc-correct.h"
#include "germline-normal-wgs/germline-normal-wgs.h"
#include "likelihood/likelihood.h"
#include "merge-segments/merge-segments-main.h"
#include "one-sample-logr/one-sample-logr.h"
#include "purity-ploidy-search/purity-ploidy-search.h"
#include "seg-to-vcf/seg-to-vcf.h"
#include "segmentation/pscbs.h"
#include "somatic-tumor-normal-wgs/somatic-tumor-normal-wgs.h"
#include "two-sample-logr/two-sample-logr.h"

namespace xoos::cnc {

namespace opt = cli_opt_name;

constexpr std::array<const char*, 11> kHiddenSubcommands = {
    kAugmentBaitsSubcommand,
    kCalculateCoverageSubcommand,
    kGCCorrectSubcommand,
    kDenoiseSubcommand,
    kTwoSampleLogRSubcommand,
    kOneSampleLogRSubcommand,
    kSegmentationSubcommand,
    kPurityPloidySearchSubcommand,
    kPredictGermlineCNVSubcommand,
    kMergeSegmentsSubcommand,
    kSegToVcfSubcommand,
};

bool IsHiddenSubcommand(const std::string& name) {
  return std::ranges::find(kHiddenSubcommands, name) != std::ranges::end(kHiddenSubcommands);
}

void ApplyHiddenGroup(CLI::App* subcommand) {
  if (subcommand != nullptr && IsHiddenSubcommand(subcommand->get_name())) {
    subcommand->group("");
  }
}

void ConfigureSomaticLikelihoodOptions(CopyNumberCallerOptions& options) {
  options.likelihood_options.mode = LikelihoodMode::kSomatic;
  options.likelihood_options.mapq_cutoff_for_calls = kLikelihoodDefaultSomaticMAPQCutoffForCalls;
  using enum segmentation::SegmentType;
  options.likelihood_options.input_segments_type = kBaf;
  options.likelihood_options.output_segments_type =
      options.bafs_fname.has_value() ? kSomaticWithBafLikelihood : kSomaticNoBafLikelihood;
}

void ConfigureGermlineLikelihoodOptions(CopyNumberCallerOptions& options) {
  options.sample_metadata_options.purity = 0.99;
  options.sample_metadata_options.ploidy = 2;
  options.likelihood_options.mode = LikelihoodMode::kGermline;
  options.likelihood_options.input_segments_type = segmentation::SegmentType::kLogROnly;
  options.likelihood_options.output_segments_type = segmentation::SegmentType::kGermlineLikelihood;
}

void ConfigureSomaticTumorNormalWGSOptions(cli::ConstAppPtr app, CopyNumberCallerOptions& options) {
  options.augment_baits_options.generate_whole_genome_baits = true;
  options.mode = CopyNumberCallerModes::kSomaticTumorNormalWGS;
  options.segmentation_options.segmentation_mode = segmentation::SegmentationMode::kSomatic;
  options.segmentation_options.undo_logr_segments = true;
  // Restore somatic-specific defaults for bindings shared with the germline subcommand.
  // CLI11 writes default_val to the binding at registration time (last writer wins), so
  // whichever subcommand was registered last clobbers the other's defaults.
  if (app->count(opt::kIntervalSize) == 0) {
    options.augment_baits_options.whole_genome_interval_size = kAugmentBaitsDefaultWholeGenomeIntervalSize;
  }
  if (app->count(opt::kMinIntervalMappability) == 0) {
    options.augment_baits_options.whole_genome_min_mappability = kAugmentBaitsDefaultWholeGenomeMinMappability;
  }
  if (app->count(opt::kMinVcfNormalDepth) == 0) {
    options.baf_filter_options.normal_sample_min_depth = kVcfParsingOptionsDefaultSomaticNormalSampleMinDepth;
  }
}

void ConfigureSomaticTumorTargetedEnrichmentOptions(CopyNumberCallerOptions& options) {
  options.mode = CopyNumberCallerModes::kSomaticTumorTargetedEnrichment;
  options.segmentation_options.segmentation_mode = segmentation::SegmentationMode::kSomatic;
}

void ConfigureGermlineNormalWGSOptions(cli::ConstAppPtr app, CopyNumberCallerOptions& options) {
  options.augment_baits_options.generate_whole_genome_baits = true;
  options.mode = CopyNumberCallerModes::kGermlineNormalWGS;
  options.sample_metadata_options.ploidy = 2;
  options.sample_metadata_options.purity = 0.99;
  options.likelihood_options.mode = LikelihoodMode::kGermline;
  options.segmentation_options.cbs_method = segmentation::CbsMaxTMethod::kFast;
  options.segmentation_options.truncate_outliers = true;
  options.segmentation_options.segmentation_mode = segmentation::SegmentationMode::kGermline;
  options.baf_filter_options.normal_sample_min_baf = kVcfParsingOptionsDefaultGermlineNormalSampleMinBAF;
  options.baf_filter_options.normal_sample_max_baf = kVcfParsingOptionsDefaultGermlineNormalSampleMaxBAF;
  // Restore germline-specific defaults for bindings shared with the somatic subcommand.
  // CLI11 writes default_val to the binding at registration time (last writer wins), so
  // whichever subcommand was registered last clobbers the other's defaults.
  if (app->count(opt::kIntervalSize) == 0) {
    options.augment_baits_options.whole_genome_interval_size = kAugmentBaitsGermlineWGSDefaultIntervalSize;
  }
  if (app->count(opt::kMinIntervalMappability) == 0) {
    options.augment_baits_options.whole_genome_min_mappability =
        kAugmentBaitsGermlineWGSDefaultWholeGenomeMinMappabiilty;
  }
  if (app->count(opt::kMinVcfNormalDepth) == 0) {
    options.baf_filter_options.normal_sample_min_depth = kVcfParsingOptionsDefaultGermlineNormalSampleMinDepth;
  }
  options.likelihood_options.mapq_cutoff_for_calls_is_user_set = app->count(opt::kMinMapqForCalls) > 0;
}

void CopyNumberCallerCliPreCallback(cli::ConstAppPtr app, CopyNumberCallerOptions& options) {
  file::CreateWritableDirectory(fs::absolute(options.output_dir));
  // update sample ID when sample ID is not exposed as an option to the user
  const auto* sample_id_opt = app->get_option_no_throw(opt::kSampleId);
  if (sample_id_opt == nullptr) {
    if (options.sample_metadata_options.normal_sample_name.has_value() &&
        options.sample_metadata_options.tumor_sample_name.has_value()) {
      // in the case of tumor-normal paired analysis, use the format normal_tumor as the sample ID
      const auto& normal_name = options.sample_metadata_options.normal_sample_name.value();
      const auto& tumor_name = options.sample_metadata_options.tumor_sample_name.value();
      if (!normal_name.empty() && !tumor_name.empty()) {
        options.sample_metadata_options.sample_id = fmt::format("{}_{}", normal_name, tumor_name);
      }
    } else if (options.sample_metadata_options.normal_sample_name.has_value()) {
      // in the case of normal-only analysis, use the normal sample name as the sample ID
      const auto& normal_name = options.sample_metadata_options.normal_sample_name.value();
      if (!normal_name.empty()) {
        options.sample_metadata_options.sample_id = normal_name;
      }
    } else {
      options.sample_metadata_options.sample_id = kSampleMetadataOptionsDefaultSampleId;
    }
  }
  const auto* reference_opt = app->get_option_no_throw(opt::kReference);
  if (reference_opt != nullptr && reference_opt->count() > 0) {
    options.reference_genome_fai_fname = options.reference_genome_fname;
    options.reference_genome_fai_fname += ".fai";
  }
  options.command_line_info = cli::GetCommandLineInfo(app);
}

static void AddHiddenSubcommands(cli::AppPtr app, CopyNumberCallerOptionsPtr& options) {
  auto* augment_baits_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kAugmentBaitsSubcommand,
      DefineAugmentBaitsSubcommandOptions,
      options,
      AugmentBaitsMain,
      "Generate or augment bait intervals.",
      [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) { CopyNumberCallerCliPreCallback(app, *opts); });
  ApplyHiddenGroup(augment_baits_sub_app);

  auto* calculate_coverage_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kCalculateCoverageSubcommand,
      DefineCalculateCoverageSubcommandOptions,
      options,
      CalculateCoverageMain,
      "Calculate coverage and mapping quality metrics.",
      [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) { CopyNumberCallerCliPreCallback(app, *opts); });
  ApplyHiddenGroup(calculate_coverage_sub_app);

  auto* gc_correct_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kGCCorrectSubcommand,
      DefineGCCorrectSubcommandOptions,
      options,
      GCCorrectMain,
      "Correct coverage for GC bias.",
      [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) { CopyNumberCallerCliPreCallback(app, *opts); });
  ApplyHiddenGroup(gc_correct_sub_app);

  auto* denoise_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kDenoiseSubcommand,
      DefineDenoiseSubcommandOptions,
      options,
      DenoiseMain,
      "Denoise coverage and compute log ratios.",
      [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) { CopyNumberCallerCliPreCallback(app, *opts); });
  ApplyHiddenGroup(denoise_sub_app);

  auto* two_sample_logr_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kTwoSampleLogRSubcommand,
      DefineTwoSampleLogRSubcommandOptions,
      options,
      TwoSampleLogRMain,
      "Compute tumor/normal log ratios.",
      [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) { CopyNumberCallerCliPreCallback(app, *opts); });
  ApplyHiddenGroup(two_sample_logr_sub_app);

  auto* one_sample_logr_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kOneSampleLogRSubcommand,
      DefineOneSampleLogRSubcommandOptions,
      options,
      OneSampleLogRMain,
      "Compute log ratios from a single sample.",
      [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) { CopyNumberCallerCliPreCallback(app, *opts); });
  ApplyHiddenGroup(one_sample_logr_sub_app);

  auto* segmentation_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kSegmentationSubcommand,
      DefineSegmentationSubcommandOptions,
      options,
      segmentation::ParentSpecificBinarySegmentationMain,
      "Segment log ratio data.",
      [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) { CopyNumberCallerCliPreCallback(app, *opts); });
  ApplyHiddenGroup(segmentation_sub_app);

  auto* purity_ploidy_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kPurityPloidySearchSubcommand,
      DefinePurityPloidySearchSubcommandOptions,
      options,
      PurityPloidySearchMain,
      "Search for purity/ploidy solutions.",
      [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) { CopyNumberCallerCliPreCallback(app, *opts); });
  ApplyHiddenGroup(purity_ploidy_sub_app);

  auto* predict_germline_sub_app =
      cli::AddSubcommand<CopyNumberCallerOptions>(app,
                                                  kPredictGermlineCNVSubcommand,
                                                  DefinePredictGermlineCNVSubcommandOptions,
                                                  options,
                                                  CalculateLikelihoodsMain,
                                                  "Predict germline CNVs from segments.",
                                                  [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) {
                                                    ConfigureGermlineLikelihoodOptions(*opts);
                                                    CopyNumberCallerCliPreCallback(app, *opts);
                                                  });
  ApplyHiddenGroup(predict_germline_sub_app);

  auto* merge_segments_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kMergeSegmentsSubcommand,
      DefineMergeSegmentsSubcommandOptions,
      options,
      MergeSegmentsMain,
      "Merge adjacent segments.",
      [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) { CopyNumberCallerCliPreCallback(app, *opts); });
  ApplyHiddenGroup(merge_segments_sub_app);

  auto* seg_to_vcf_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kSegToVcfSubcommand,
      DefineSegToVcfSubcommandOptions,
      options,
      segmentation::WriteSegmentsToVcfMain,
      "Convert SEG files to VCF.",
      [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) { CopyNumberCallerCliPreCallback(app, *opts); });
  ApplyHiddenGroup(seg_to_vcf_sub_app);
}

void AddSubcommands(cli::AppPtr app, CopyNumberCallerOptionsPtr& options) {
  AddHiddenSubcommands(app, options);

  cli::AddSubcommand<CopyNumberCallerOptions>(app,
                                              kPredictSomaticCNASubcommand,
                                              DefinePredictSomaticCNASubcommandOptions,
                                              options,
                                              CalculateLikelihoodsMain,
                                              "Predict somatic CNAs from segments.",
                                              [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) {
                                                ConfigureSomaticLikelihoodOptions(*opts);
                                                CopyNumberCallerCliPreCallback(app, *opts);
                                              });

  cli::AddSubcommand<CopyNumberCallerOptions>(app,
                                              kSomaticTumorNormalWGSSubcommand,
                                              DefineSomaticTumorNormalWGSOptions,
                                              options,
                                              SomaticTumorNormalWGSMain,
                                              "Run somatic tumor/normal WGS pipeline.",
                                              [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) {
                                                ConfigureSomaticTumorNormalWGSOptions(app, *opts);
                                                CopyNumberCallerCliPreCallback(app, *opts);
                                              });

  cli::AddSubcommand<CopyNumberCallerOptions>(app,
                                              kGermlineNormalWGSSubcommand,
                                              DefineGermlineNormalWGSOptions,
                                              options,
                                              GermlineNormalWGSMain,
                                              "Run germline normal WGS pipeline.",
                                              [](cli::ConstAppPtr app, const CopyNumberCallerOptionsPtr& opts) {
                                                ConfigureGermlineNormalWGSOptions(app, *opts);
                                                CopyNumberCallerCliPreCallback(app, *opts);
                                              });
}

s32 CopyNumberCallerCliMain(s32 argc, char** argv) {
  auto app = cli::SetupDefaultCli(PROGRAM_NAME, VERSION);
  app->require_subcommand(1);
  auto options = std::make_shared<CopyNumberCallerOptions>();
  AddSubcommands(app.get(), options);
  return xoos::cli::RunCli(app.get(), argc, argv);
}
}  // namespace xoos::cnc
