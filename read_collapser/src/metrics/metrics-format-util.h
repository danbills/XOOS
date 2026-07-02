#pragma once

#include <string>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::read_collapser {

// kNA represents metric fields that are not applicable (e.g. consensus metrics will be NA for duplicate marking mode).
constexpr std::string kNA = "NA";

// Returns the percentage of value/total as a formatted string with the given precision.
std::string CalculatePercentage(u64 value, u64 total, u8 precision = 2);

// Returns a formatted TSV row with metric name, value, percentage, and denominator name.
// If metric_not_na is false, value and percentage fields are set to kNA.
vec<std::string> FormatRow(const std::string& metric_name,
                           u64 value,
                           u64 denominator,
                           const std::string& denominator_name,
                           bool metric_not_na = true);

}  // namespace xoos::read_collapser
