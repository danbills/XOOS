#include "model-resolver.h"

#include <fmt/format.h>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include "core/cli-option-names.h"

namespace xoos::svc {

using enum AlignerType;

// Returns the formatted aligner name, rejecting kCustom since it requires explicit model paths.
static std::string AlignerToModelString(const AlignerType aligner) {
  if (aligner == kCustom) {
    throw error::Error(fmt::format("Cannot resolve model for {} aligner — provide {} explicitly",
                                   enum_util::FormatEnumName(aligner),
                                   cli_opt_name::kModel));
  }
  return enum_util::FormatEnumName(aligner);
}

fs::path ResolveTumorNormalModel(const SampleType sample_type, const AlignerType aligner) {
  // No pre-trained giraffe models exist for tumor-normal-wgs — reject rather than silently
  // resolving a mismatched model.
  if (aligner == kGiraffe) {
    throw error::Error(
        fmt::format("Aligner '{}' is not supported for the tumor-normal-wgs workflow because no pre-trained {} model "
                    "is available. Use '{} {}', or '{} {}' with an explicit {}.",
                    enum_util::FormatEnumName(kGiraffe),
                    enum_util::FormatEnumName(kGiraffe),
                    cli_opt_name::kAligner,
                    enum_util::FormatEnumName(kBwa),
                    cli_opt_name::kAligner,
                    enum_util::FormatEnumName(kCustom),
                    cli_opt_name::kModel));
  }
  const auto sample_type_str = enum_util::FormatEnumName(sample_type);
  const auto path =
      fmt::format("/resources/model-tumor-normal-wgs-{}-{}.txt.gz", sample_type_str, AlignerToModelString(aligner));

  Logging::Info("Resolved model: {} based on aligner {} and sample_type {}",
                path,
                enum_util::FormatEnumName(aligner),
                sample_type_str);
  return fs::path{path};
}

// Shared logic for resolving germline SNV/indel model paths.
// Validates the aligner and constructs paths from the templates.
// Each template must contain two {} placeholders: run_type and aligner.
static std::pair<fs::path, fs::path> ResolveGermlineModelPair(const RunType run_type,
                                                              const AlignerType aligner,
                                                              const std::string_view snv_template,
                                                              const std::string_view indel_template) {
  if (aligner == kCustom) {
    throw error::Error(fmt::format("Cannot resolve germline models for {} aligner — provide {} and {} explicitly",
                                   enum_util::FormatEnumName(aligner),
                                   cli_opt_name::kSnvModel,
                                   cli_opt_name::kIndelModel));
  }

  const auto run_type_str = enum_util::FormatEnumName(run_type);
  const auto aligner_str = AlignerToModelString(aligner);

  const auto snv_path = fmt::format(fmt::runtime(snv_template), run_type_str, aligner_str);
  const auto indel_path = fmt::format(fmt::runtime(indel_template), run_type_str, aligner_str);

  Logging::Info("Resolved SNV model: {} for aligner {} and run type {}", snv_path, aligner_str, run_type_str);
  Logging::Info("Resolved indel model: {} for aligner {} and run type {}", indel_path, aligner_str, run_type_str);

  return {fs::path{snv_path}, fs::path{indel_path}};
}

std::pair<fs::path, fs::path> ResolveGermlineModels(const RunType run_type, const AlignerType aligner) {
  return ResolveGermlineModelPair(
      run_type, aligner, "/resources/model-germline-{}-{}-snv.txt.gz", "/resources/model-germline-{}-{}-indel.txt.gz");
}

std::pair<fs::path, fs::path> ResolveGermlineMultiSampleModels(const RunType run_type, const AlignerType aligner) {
  return ResolveGermlineModelPair(run_type,
                                  aligner,
                                  "/resources/model-germline-{}-{}-multisample-snv.txt.gz",
                                  "/resources/model-germline-{}-{}-multisample-indel.txt.gz");
}

std::string TumorNormalModelThresholdsKey(const SampleType sample_type, const AlignerType aligner) {
  return fmt::format("{}-{}", enum_util::FormatEnumName(sample_type), enum_util::FormatEnumName(aligner));
}

std::optional<ModelThresholds> ResolveTumorNormalThresholds(const SampleType sample_type,
                                                            const AlignerType aligner,
                                                            const StrMap<ModelThresholds>& model_thresholds) {
  if (aligner == kCustom) {
    throw error::Error(fmt::format("Cannot resolve thresholds for {} aligner — provide explicit score thresholds",
                                   enum_util::FormatEnumName(aligner)));
  }

  if (model_thresholds.empty()) {
    return std::nullopt;
  }

  // Giraffe is rejected by ResolveTumorNormalModel before this is reached, so only bwa keys are resolved here.
  const auto key = TumorNormalModelThresholdsKey(sample_type, aligner);
  const auto it = model_thresholds.find(key);
  if (it == model_thresholds.end()) {
    Logging::Warn("No model_thresholds entry for '{}', using top-level config defaults", key);
    return std::nullopt;
  }

  Logging::Info("Resolved scoring thresholds for '{}': snv_min_ml_score={}, indel_min_ml_score={}",
                key,
                it->second.snv_min_ml_score,
                it->second.indel_min_ml_score);
  return it->second;
}

}  // namespace xoos::svc
