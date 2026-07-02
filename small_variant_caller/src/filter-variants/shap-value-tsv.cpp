#include "shap-value-tsv.h"

#include <fmt/format.h>

#include "xoos/error/error.h"

namespace xoos::svc {

void WriteShapValueTsvHeader(LockedTsvWriter& writer,
                             const vec<std::string>& feature_names,
                             const std::optional<io::CommandLineInfo>& cmd_info) {
  if (feature_names.empty()) {
    throw error::Error("Feature names cannot be empty when writing SHAP value TSV header.");
  }
  if (cmd_info.has_value()) {
    writer.AppendMetadata(cmd_info.value());
  }
  vec<std::string> header = {"VARIANT", "GENOTYPE", "PRED_ML"};
  header.insert(header.end(), feature_names.begin(), feature_names.end());
  header.emplace_back("BASE_VALUE");
  writer.AppendRow(header);
}

vec<std::string> AssembleShapValueTsvRow(const VariantId& vid, const PredictionScore& score) {
  if (score.shap_values.empty()) {
    throw error::Error("SHAP values cannot be empty when writing SHAP value TSV row for variant {}", vid.ToString());
  }
  vec<std::string> row = {vid.ToString(), GenotypeToString(score.genotype), fmt::format("{:.17g}", score.probability)};
  for (const auto& shap_value : score.shap_values) {
    row.emplace_back(fmt::format("{:.17g}", shap_value));
  }
  return row;
}

}  // namespace xoos::svc
