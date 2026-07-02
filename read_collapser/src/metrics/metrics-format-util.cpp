#include "metrics/metrics-format-util.h"

#include <fmt/format.h>

namespace xoos::read_collapser {

std::string CalculatePercentage(const u64 value, const u64 total, const u8 precision) {
  if (total == 0) {
    return kNA;
  }
  const f64 percentage = static_cast<f64>(value) / static_cast<f64>(total) * 100.0;
  return fmt::format("{:.{}f}", percentage, precision);
}

vec<std::string> FormatRow(const std::string& metric_name,
                           const u64 value,
                           const u64 denominator,
                           const std::string& denominator_name,
                           const bool metric_not_na) {
  if (metric_not_na) {
    vec<std::string> output{metric_name, std::to_string(value)};
    if (denominator_name != kNA && denominator != 0) {
      output.emplace_back(CalculatePercentage(value, denominator));
    } else {
      output.push_back(kNA);
    }
    output.push_back(denominator_name);
    return output;
  }
  return {metric_name, kNA, kNA, denominator_name};
}

}  // namespace xoos::read_collapser
