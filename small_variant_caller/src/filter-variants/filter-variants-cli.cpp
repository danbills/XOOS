#include "filter-variants-cli.h"

#include <algorithm>
#include <array>
#include <string>

#include <fmt/format.h>

#include <xoos/cli/cli.h>
#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/file-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/error/error.h>

#include "../core/workflow.h"
#include "../util/model-resolver.h"
#include "compute-bam-features/compute-bam-features-cli.h"
#include "compute-vcf-features/compute-vcf-features-cli.h"
#include "core/cli-option-names.h"
#include "util/cli-util.h"
#include "util/model-resolver.h"

namespace xoos::svc {

using cli::AddOptionalEnumOption;
using cli::AddThreadCountOption;
using enum Workflow;

// "custom" workflow not supported in filter_variants
constexpr std::array kSupportedWorkflows = {
    kGermline, kGermlineMultiSample, kTumorOnlyTe, kTumorNormalWgs, kGermlineTagging};

// CLI option group names
constexpr auto* kCliOptGroupBamFeatureExtraction{"BAM feature extraction"};
constexpr auto* kCliOptGroupVcfFeatureExtraction{"VCF feature extraction"};

// Default values for CLI options
constexpr u32 kDefaultMaxVcfRegionSizePerThread = 64'000;
constexpr auto* kDefaultSdChrName = "chr1";
constexpr auto* kDefaultVcfOutput = "output.vcf.gz";
// Default pre-trained model path(s) for each workflow
constexpr auto* kDefaultTumorOnlyTeModelPath = "/resources/model-ffpe-bwa.txt.gz";

/**
 * @brief Helper function to add core CLI options for subcommands. CLI options added here are common across all
 * subcommands.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
static void AddCommonOptions(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  AddWarnAsErrorOption(sub);

  AddThreadCountOption(sub, cli_opt_name::kThreads, params->threads);

  sub->add_option(cli_opt_name::kConfig, params->config_file, "Path to config JSON file")->check(CLI::ExistingFile);

  auto* const bam_opt = sub->add_option(cli_opt_name::kBamInput,
                                        params->bam_files,
                                        "Path(s) to input BAM file(s), produced by GATK HaplotypeCaller/Mutect2")
                            ->required();
  CheckIndexedBamFile(bam_opt);

  auto* const vcf_opt = sub->add_option(cli_opt_name::kVcfInput,
                                        params->vcf_file,
                                        "Path to input VCF file, produced by GATK Mutect2/HaplotypeCaller")
                            ->required();
  CheckIndexedVcfFile(vcf_opt);

  auto* const genome_opt =
      sub->add_option(cli_opt_name::kGenome, params->genome, "Path to indexed FASTA file for reference genome")
          ->required();
  CheckIndexedFastaFile(genome_opt);

  auto* const regions_opt =
      sub->add_option(cli_opt_name::kTargetRegions, params->bed_file, "Path to BED file for 0-based target regions");
  CheckBedFile(regions_opt);

  sub->add_option(cli_opt_name::kOutputDir, params->output_dir, "Path to output directory")->default_val(".");

  const auto vcf_output_desc = fmt::format("Path to output VCF file, relative to {}", cli_opt_name::kOutputDir);
  sub->add_option(cli_opt_name::kVcfOutput, params->vcf_output, vcf_output_desc)
      ->default_val(kDefaultVcfOutput)
      ->transform([&params](const std::string& path) { return (params->output_dir / fs::path(path)).string(); })
      ->check(CLI::NonexistentPath);

  sub->add_option(cli_opt_name::kMaxBamRegionSizePerThread,
                  params->max_bam_region_size_per_thread,
                  "Maximum BAM region size per thread during feature extraction (inclusive)")
      ->default_val(kDefaultMaxBamRegionSizePerThread)
      ->check(CLI::PositiveNumber);

  sub->add_option(cli_opt_name::kMaxVcfRegionSizePerThread,
                  params->max_vcf_region_size_per_thread,
                  "Maximum VCF region size per thread during variant filtration (inclusive)")
      ->default_val(kDefaultMaxVcfRegionSizePerThread)
      ->check(CLI::PositiveNumber);

  sub->add_option(cli_opt_name::kSdChrName, params->sd_chr_name, "autosome name for sex determination")
      ->default_val(kDefaultSdChrName);

  auto* const par_bed_x_opt = sub->add_option(
      cli_opt_name::kParBedX, params->par_bed_x, "Path to chromosome X pseudoautosommal region BED file");
  CheckBedFile(par_bed_x_opt);

  auto* const par_bed_y_opt = sub->add_option(
      cli_opt_name::kParBedY, params->par_bed_y, "Path to chromosome Y pseudoautosommal region BED file");
  CheckBedFile(par_bed_y_opt);
}

/**
 * @brief Helper function to add CLI options for outputting SHAP values.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
static void AddOutputShapValueTsvOption(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  sub->add_option(cli_opt_name::kOutputShapValueTsv,
                  params->shap_value_tsv,
                  "Path to the SHAP value output TSV file, relative to output directory. Intended to be applied to "
                  "target regions due to slow runtime.")
      ->needs(cli_opt_name::kTargetRegions);
}

/**
 * @brief Helper function to add CLI options for outputting SHAP values for SNVs and indels separately. This is intended
 * to be used for germline workflow where separate SHAP value output for SNVs and indels is applicable.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
static void AddOutputSnvIndelShapValueTsvOptions(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  sub->add_option(cli_opt_name::kOutputSnvShapValueTsv,
                  params->snv_shap_value_tsv,
                  "Path to the SHAP value output TSV file for SNVs, relative to output directory. Intended to be "
                  "applied to target regions due to slow runtime.")
      ->needs(cli_opt_name::kTargetRegions);

  sub->add_option(cli_opt_name::kOutputIndelShapValueTsv,
                  params->indel_shap_value_tsv,
                  "Path to the SHAP value output TSV file for indels, relative to output directory. Intended to be "
                  "applied to target regions due to slow runtime.")
      ->needs(cli_opt_name::kTargetRegions);
}

/**
 * @brief Add an --aligner option that accepts only the aligners supported by tumor-normal-wgs.
 * @details tumor-normal-wgs only ships pre-trained BWA models, so 'bwa' (auto-select) and 'custom'
 *          (explicit --model) are the only valid choices. 'giraffe' is intentionally not accepted —
 *          it is rejected at parse time and is not advertised in the help text — to avoid the
 *          confusing combination of '--aligner giraffe' with an explicit '--model'.
 * @param sub Subcommand application pointer where the option is added.
 * @param params CLI parameters shared pointer to store the parsed aligner value.
 */
static void AddTumorNormalWgsAlignerOption(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  static constexpr std::array kSupportedAligners = {AlignerType::kBwa, AlignerType::kCustom};
  const auto supported = fmt::format(
      "{{{}, {}}}", enum_util::FormatEnumName(AlignerType::kBwa), enum_util::FormatEnumName(AlignerType::kCustom));
  sub->add_option_function<std::string>(
         cli_opt_name::kAligner,
         [params, supported](const std::string& input) {
           const auto parsed = enum_util::ParseEnumName<AlignerType>(input);
           if (!parsed || std::ranges::find(kSupportedAligners, *parsed) == kSupportedAligners.end()) {
             throw CLI::ValidationError(
                 fmt::format("{}: Check {} value in {} FAILED. 'giraffe' is not supported for "
                             "the tumor-normal-wgs workflow because no pre-trained giraffe model "
                             "is available; use 'custom' with an explicit {} instead.",
                             cli_opt_name::kAligner,
                             input,
                             supported,
                             cli_opt_name::kModel));
           }
           params->aligner = *parsed;
         },
         fmt::format("Aligner used for BAM alignment. Set to 'bwa' to auto-select the model based on --sample-type, "
                     "or leave as 'custom' and provide --model explicitly. Accepted values: {}",
                     supported))
      ->run_callback_for_default()
      ->default_val(enum_util::FormatEnumName(AlignerType::kCustom));
}

/**
 * @brief Helper function to add tumor-normal-wgs specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddTumorNormalWgsSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kTumorNormalWgs));
  const auto defaults = SVCConfig(kTumorNormalWgs);

  // Model is required by default (aligner defaults to 'custom').
  // When --aligner is explicitly set to 'bwa', the model is auto-resolved.
  auto* const model_opt = sub->add_option(cli_opt_name::kModel, params->model)
                              ->description(
                                  "Path to input model file. Required by default. "
                                  "When --aligner is set to 'bwa', the model is automatically "
                                  "selected based on --sample-type");
  CheckNonEmptyFile(model_opt);

  cli::AddEnumOption(sub,
                     cli_opt_name::kSampleType,
                     params->sample_type,
                     "Sample preparation type for model selection (used when --aligner is 'bwa')",
                     SampleType::kFfpe);

  AddTumorNormalWgsAlignerOption(sub, params);

  sub->add_option(cli_opt_name::kTumorSampleName, params->tumor_sample_name, "tumor sample name for read groups")
      ->required();

  sub->add_option(
         cli_opt_name::kSnvMinMlScore, params->snv_min_ml_score, "Minimum ML score threshold (inclusive) for SNVs")
      ->default_val(defaults.snv_min_ml_score)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kIndelMinMlScore,
                  params->indel_min_ml_score,
                  "Minimum ML score threshold (inclusive) for indels")
      ->default_val(defaults.indel_min_ml_score)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMinTumorSupport,
                  params->min_tumor_support,
                  "Minimum support threshold (inclusive) for variants in the tumor sample")
      ->default_val(defaults.min_tumor_support)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMaxNormalSupport,
                  params->max_normal_support,
                  "Maximum support threshold (inclusive) for variants in the normal sample")
      ->default_val(defaults.max_normal_support)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMinTumorAf,
                  params->min_tumor_af,
                  "Minimum allele frequency threshold (inclusive) for variants in the tumor sample")
      ->default_val(defaults.min_tumor_af)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMinDpRatio,
                  params->min_dp_ratio,
                  "Minimum ratio (inclusive) of DP to chromosome's median DP for variants in the tumor sample")
      ->default_val(defaults.min_dp_ratio)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMaxIndelSize, params->max_indel_size, "Maximum indel size in base pairs (inclusive)")
      ->default_val(defaults.max_indel_size)
      ->check(CLI::PositiveNumber);

  AddOutputShapValueTsvOption(sub, params);
}

/**
 * @brief Helper function to add tumor-only-te specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddTumorOnlyTeSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kTumorOnlyTe));
  const auto defaults = SVCConfig(kTumorOnlyTe);

  auto* const model_opt = sub->add_option(cli_opt_name::kModel, params->model)
                              ->description("Path to input model file")
                              ->default_val(kDefaultTumorOnlyTeModelPath);
  CheckNonEmptyFile(model_opt);

  auto* const blocklist_opt = sub->add_option(cli_opt_name::kBlocklist,
                                              params->block_list,
                                              "Text file listing variants to skip. Variants are represented as "
                                              "`chr_pos_ref_alt`, one per line. Variant position is 1-based.");
  CheckNonEmptyFile(blocklist_opt);

  auto* const hotspot_opt = sub->add_option(
      cli_opt_name::kHotspotVcf, params->hotspot_list, "A VCF file containing a list of hotspot variants");
  CheckNonEmptyFile(hotspot_opt);

  auto* const forcecall_opt =
      sub->add_option(cli_opt_name::kForcecallBed,
                      params->forcecall_list,
                      "A 0-based BED file with forced call positions to be included in the output");
  CheckBedFile(forcecall_opt);

  sub->add_option(cli_opt_name::kMinAltCounts,
                  params->min_alt_counts,
                  "Minimum number of alt counts for retaining a variant (inclusive)")
      ->default_val(defaults.min_alt_counts)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMinAf,
                  params->min_allele_freq_threshold,
                  "Minimum allowed variant allele frequency (inclusive)")
      ->default_val(defaults.min_af)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMinPhasedAf,
                  params->min_phased_allele_freq,
                  "Minimum allowed phased allele frequency for a variant (inclusive)")
      ->default_val(defaults.min_phased_af)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMaxPhasedAf,
                  params->max_phased_allele_freq,
                  "Maximum allowed phased allele frequency for a variant (inclusive)")
      ->default_val(defaults.max_phased_af)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMinWeightedCounts,
                  params->weighted_counts_threshold,
                  "Threshold for weighted counts (inclusive). Variants must score at least this much to be included")
      ->default_val(defaults.min_weighted_counts)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(
         cli_opt_name::kHotspotMinWeightedCounts,
         params->hotspot_weighted_counts_threshold,
         "Threshold for hotspot weighted counts (inclusive). Variants must score at least this much to be included")
      ->default_val(defaults.hotspot_min_weighted_counts)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMinMlScore,
                  params->ml_threshold,
                  "Threshold for ML score (inclusive). Variants must score at least this much to be included")
      ->default_val(defaults.min_ml_score)
      ->check(kCliRangeFraction);

  sub->add_option(
         cli_opt_name::kHotspotMinMlScore,
         params->hotspot_ml_threshold,
         "Threshold for ML score in hotspots (inclusive). Variants must score at least this much to be included")
      ->default_val(defaults.hotspot_min_ml_score)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kPhased, params->phased, "Call phased variants in VCF")->default_val(defaults.phased);

  AddOutputShapValueTsvOption(sub, params);
}

/**
 * @brief Helper function to add germline-tagging specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddGermlineTaggingSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging));

  // no default model for germline-tagging, user must provide one
  auto* const opt =
      sub->add_option(cli_opt_name::kModel, params->model)->description("Path to input model file")->required();
  CheckNonEmptyFile(opt);

  AddOutputShapValueTsvOption(sub, params);
}

/**
 * @brief Helper function to add germline-common CLI options (--aligner, --run-type, --snv-model, --indel-model).
 *        Used by both germline and germline-multi-sample subcommands.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
static void AddGermlineModelSelectionOptions(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  cli::AddEnumOption(sub,
                     cli_opt_name::kAligner,
                     params->aligner,
                     "Aligner used for BAM alignment. Set to 'bwa' or 'giraffe' to auto-select models, "
                     "or leave as 'custom' and provide --snv-model and --indel-model explicitly",
                     AlignerType::kCustom);

  cli::AddEnumOption(
      sub, cli_opt_name::kRunType, params->run_type, "Run type for model selection (sbxd or sbxfast)", RunType::kSbxd);

  auto* const snv_model_opt = sub->add_option(cli_opt_name::kSnvModel, params->snv_model)
                                  ->description("Path to input SNV model file. Required when --aligner is 'custom'");
  CheckNonEmptyFile(snv_model_opt);

  auto* const indel_model_opt =
      sub->add_option(cli_opt_name::kIndelModel, params->indel_model)
          ->description("Path to input indel model file. Required when --aligner is 'custom'");
  CheckNonEmptyFile(indel_model_opt);
}

/**
 * @brief Helper function to add germline specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddGermlineSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kGermline));

  AddGermlineModelSelectionOptions(sub, params);
  AddOutputSnvIndelShapValueTsvOptions(sub, params);
}

/**
 * @brief Helper function to add germline-multi-sample specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddGermlineMultiSampleSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kGermlineMultiSample));

  AddGermlineModelSelectionOptions(sub, params);
  AddOutputSnvIndelShapValueTsvOptions(sub, params);
}

void filter_variants::DefineOptions(CLI::App* const app, FilterVariantsParamPtr& params) {
  // Add a subcommand for each workflow
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    const std::string desc = fmt::format("Filter variants for the {} workflow", name);
    CLI::App* const sub = app->add_subcommand(name, desc)->fallthrough();
    // Do not apply force_callback() to subcommand options to avoid overwriting params set by other subcommands
    AddCommonOptions(sub, params);

    // Add options for VCF feature extraction
    for (auto* const opt : compute_vcf_features::AddSharedOptions(sub, params)) {
      opt->group(kCliOptGroupVcfFeatureExtraction);
    }

    const auto defaults = SVCConfig(workflow);

    // Add options for BAM feature extraction
    for (auto* const opt : compute_bam_features::AddSharedOptions(sub, params, defaults)) {
      opt->group(kCliOptGroupBamFeatureExtraction);
    }

    cli::AddEnumOption(sub,
                       cli_opt_name::kNormalizeFeatures,
                       params->normalize_features,
                       "Normalize read-depth related feature values",
                       defaults.normalize_features);
  }
  app->require_subcommand(kMinSubcommands, kMaxSubcommands);

  // Add workflow-specific CLI options for each subcommand
  AddTumorNormalWgsSpecificOptions(app, params);
  AddTumorOnlyTeSpecificOptions(app, params);
  AddGermlineTaggingSpecificOptions(app, params);
  AddGermlineSpecificOptions(app, params);
  AddGermlineMultiSampleSpecificOptions(app, params);

  // hide the `germline-tagging` subcommand by assigning it to an empty string group
  app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging))->group("");
}

/**
 * @brief Helper function to apply workflow config to tumor-normal-wgs unique CLI parameters.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
static void ApplyConfigToTumorNormalWgsUniqueOptions(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  const bool is_not_tn_wgs = (params->config.workflow != kTumorNormalWgs);
  if (is_not_tn_wgs || sub->count(cli_opt_name::kSnvMinMlScore) == 0) {
    params->snv_min_ml_score = params->config.snv_min_ml_score;
  }
  if (is_not_tn_wgs || sub->count(cli_opt_name::kIndelMinMlScore) == 0) {
    params->indel_min_ml_score = params->config.indel_min_ml_score;
  }
  if (is_not_tn_wgs || sub->count(cli_opt_name::kMinTumorSupport) == 0) {
    params->min_tumor_support = params->config.min_tumor_support;
  }
  if (is_not_tn_wgs || sub->count(cli_opt_name::kMaxNormalSupport) == 0) {
    params->max_normal_support = params->config.max_normal_support;
  }
  if (is_not_tn_wgs || sub->count(cli_opt_name::kMinTumorAf) == 0) {
    params->min_tumor_af = params->config.min_tumor_af;
  }
  if (is_not_tn_wgs || sub->count(cli_opt_name::kMaxIndelSize) == 0) {
    params->max_indel_size = params->config.max_indel_size;
  }
  if (is_not_tn_wgs || sub->count(cli_opt_name::kMinDpRatio) == 0) {
    params->min_dp_ratio = params->config.min_dp_ratio;
  }
}

/**
 * @brief Helper function to apply workflow config to tumor-only-te unique CLI parameters.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
static void ApplyConfigToTumorOnlyTeUniqueOptions(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  const bool is_not_to_te = (params->config.workflow != kTumorOnlyTe);
  if (is_not_to_te || sub->count(cli_opt_name::kMinAltCounts) == 0) {
    params->min_alt_counts = params->config.min_alt_counts;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kMinAf) == 0) {
    params->min_allele_freq_threshold = params->config.min_af;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kMinPhasedAf) == 0) {
    params->min_phased_allele_freq = params->config.min_phased_af;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kMaxPhasedAf) == 0) {
    params->max_phased_allele_freq = params->config.max_phased_af;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kMinWeightedCounts) == 0) {
    params->weighted_counts_threshold = params->config.min_weighted_counts;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kHotspotMinWeightedCounts) == 0) {
    params->hotspot_weighted_counts_threshold = params->config.hotspot_min_weighted_counts;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kMinMlScore) == 0) {
    params->ml_threshold = params->config.min_ml_score;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kHotspotMinMlScore) == 0) {
    params->hotspot_ml_threshold = params->config.hotspot_min_ml_score;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kPhased) == 0) {
    params->phased = params->config.phased;
  }
}

static void SetNonGermlineDefaultModelPath(const CLI::App* const sub,
                                           const FilterVariantsParamPtr& params,
                                           const char* default_path) {
  if (sub->count(cli_opt_name::kModel) == 0) {
    params->model = default_path;
  }
}

// Resolve the tumor-normal-wgs model path from --sample-type and --aligner when --model is not
// provided. The default aligner is 'custom', which requires --model to be explicitly provided.
// When --aligner is set to 'bwa' or 'giraffe', the model is auto-resolved.
//
// @param sub The parsed CLI subcommand (used to check if --model was provided).
// @param params Shared pointer to CLI parameters containing aligner, sample_type, and model.
static void ResolveTumorNormalWgsModelPath(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  const bool model_provided = sub->count(cli_opt_name::kModel) > 0;

  if (params->aligner == AlignerType::kCustom && !model_provided) {
    throw error::Error(fmt::format("{} is required when {} is '{}'",
                                   cli_opt_name::kModel,
                                   cli_opt_name::kAligner,
                                   enum_util::FormatEnumName(AlignerType::kCustom)));
  }

  if (model_provided) {
    return;
  }

  params->model = ResolveTumorNormalModel(params->sample_type, params->aligner);
}

// Apply per-model ML score thresholds for tumor-normal-wgs after the top-level config defaults
// have been applied. Must run after ApplyConfigToTumorNormalWgsUniqueOptions so that model-specific
// thresholds take precedence over the top-level config values. Explicit CLI values
// (--snv-min-ml-score / --indel-min-ml-score) still take highest precedence.
//
// @param sub The parsed CLI subcommand (used to check if threshold options were explicitly set).
// @param params Shared pointer to CLI parameters containing aligner, sample_type, and thresholds.
static void ApplyModelThresholds(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  if (params->workflow != kTumorNormalWgs || params->aligner == AlignerType::kCustom) {
    return;
  }

  // Only apply when the model was auto-resolved (not explicitly provided via --model)
  if (sub->count(cli_opt_name::kModel) > 0) {
    return;
  }

  const auto thresholds =
      ResolveTumorNormalThresholds(params->sample_type, params->aligner, params->config.model_thresholds);
  if (thresholds.has_value()) {
    if (sub->count(cli_opt_name::kSnvMinMlScore) == 0) {
      params->snv_min_ml_score = thresholds->snv_min_ml_score;
    }
    if (sub->count(cli_opt_name::kIndelMinMlScore) == 0) {
      params->indel_min_ml_score = thresholds->indel_min_ml_score;
    }
  }
}

/**
 * @brief Resolve germline SNV/indel model paths from --run-type and --aligner when models are not provided.
 * @details The default aligner is 'custom', which requires --snv-model and --indel-model to be explicitly provided.
 *          When --aligner is set to 'bwa' or 'giraffe', models are auto-resolved using --run-type.
 */
static void ResolveGermlineModelPaths(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  const bool snv_model_provided = sub->count(cli_opt_name::kSnvModel) > 0;
  const bool indel_model_provided = sub->count(cli_opt_name::kIndelModel) > 0;

  // If both models are explicitly provided, use them regardless of aligner/run-type
  if (snv_model_provided && indel_model_provided) {
    return;
  }

  if (params->aligner == AlignerType::kCustom) {
    throw error::Error(fmt::format("{} and {} are required when {} is '{}'",
                                   cli_opt_name::kSnvModel,
                                   cli_opt_name::kIndelModel,
                                   cli_opt_name::kAligner,
                                   enum_util::FormatEnumName(AlignerType::kCustom)));
  }

  // Auto-resolve both models — partial override is not supported
  if (params->workflow == kGermline) {
    const auto [snv_path, indel_path] = ResolveGermlineModels(params->run_type, params->aligner);
    params->snv_model = snv_path;
    params->indel_model = indel_path;
  } else if (params->workflow == kGermlineMultiSample) {
    const auto [snv_path, indel_path] = ResolveGermlineMultiSampleModels(params->run_type, params->aligner);
    params->snv_model = snv_path;
    params->indel_model = indel_path;
  } else {
    throw error::Error("Unsupported workflow for germline model resolution");
  }
}

/**
 * @brief Helper function to set default model paths based on workflow if not provided via CLI.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
static void SetDefaultModelPaths(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  switch (params->workflow) {
    case kGermline:
    case kGermlineMultiSample:
      ResolveGermlineModelPaths(sub, params);
      break;
    case kTumorOnlyTe:
      SetNonGermlineDefaultModelPath(sub, params, kDefaultTumorOnlyTeModelPath);
      break;
    case kTumorNormalWgs:
      ResolveTumorNormalWgsModelPath(sub, params);
      break;
    default:
      SetNonGermlineDefaultModelPath(sub, params, "");
      break;
  }
}

/**
 * @brief Helper function to apply CLI defaults and workflow config presets to CLI parameters if the corresponding CLI
 * options were not specified.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
static void ApplyConfig(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  if (sub->count(cli_opt_name::kVcfOutput) == 0) {
    // update VCF output path to be relative to output directory
    params->vcf_output = fs::path(params->output_dir) / params->vcf_output;
    if (fs::exists(params->vcf_output)) {
      throw CLI::ValidationError(fmt::format("Output VCF file '{}' already exists", params->vcf_output.string()));
    }
  }

  // resolve SHAP value TSV paths relative to output directory (order-independent)
  const auto resolve_shap_path = [&params](std::optional<fs::path>& path) {
    if (path.has_value()) {
      path = params->output_dir / path.value();
      if (fs::exists(path.value())) {
        throw CLI::ValidationError(fmt::format("SHAP value output file '{}' already exists", path->string()));
      }
    }
  };
  resolve_shap_path(params->shap_value_tsv);
  resolve_shap_path(params->snv_shap_value_tsv);
  resolve_shap_path(params->indel_shap_value_tsv);

  // set default model paths if not provided via CLI
  SetDefaultModelPaths(sub, params);

  // apply config to compute-bam-features options
  compute_bam_features::ApplyConfig(sub, params);

  if (sub->count(cli_opt_name::kNormalizeFeatures) == 0) {
    params->normalize_features = params->config.normalize_features;
  }

  // apply config to subcommand unique options
  ApplyConfigToTumorNormalWgsUniqueOptions(sub, params);
  ApplyConfigToTumorOnlyTeUniqueOptions(sub, params);

  // apply per-model thresholds after top-level config defaults so they take precedence
  ApplyModelThresholds(sub, params);
}

/**
 * @brief Preprocess the CLI options and validate the parameters.
 * @param app CLI application pointer
 * @param params CLI parameters to be validated
 * @throws CLI::ValidationError if the parameters are invalid
 */
void filter_variants::PreCallback(const cli::ConstAppPtr app, const FilterVariantsParamPtr& params) {
  params->command_line = cli::GetCommandLineInfo(app);

  // Check which subcommand was used, set workflow and config accordingly, and apply config defaults as needed
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    if (app->got_subcommand(name)) {
      params->workflow = workflow;
      params->config = JsonToConfig(params->config_file.value_or(fs::path{}), params->workflow);
      ApplyConfig(app->get_subcommand(name), params);
      break;
    }
  }

  if (params->decode_yc != YcDecodeMethod::kNone && !IsDuplexProtocol(params->sequencing_protocol)) {
    throw CLI::ValidationError("Duplex sequencing protocol must be used when decoding YC tags");
  }
}

}  // namespace xoos::svc
