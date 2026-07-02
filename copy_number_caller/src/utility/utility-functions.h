#pragma once

#include <cmath>
#include <concepts>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include <xoos/types/fs.h>

namespace xoos::cnc {

/**
 * Checks if a string matches a "true" value (1, "TRUE", "true") or a "false" value (0, "FALSE", "false"). Throws an
 * exception if the string does not match any of these values.
 */
bool IsTrueString(const std::string& s);
std::tuple<std::string, size_t, size_t> ParseRegionString(const std::string& r);

/// Check if the contig in a region string exactly matches the given name, without allocating.
bool IsEqualContig(std::string_view region, std::string_view contig);

size_t StringToNonNegativeUIntOrThrow(const std::string& s);

std::string FloatAsString(std::floating_point auto val, int precision = 6) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) << static_cast<double>(val);
  return ss.str();
}

std::string FloatAsIntString(std::floating_point auto val) {
  return FloatAsString(val, 0);
}

std::vector<fs::path> PathsFromFile(const fs::path& fname);
}  // namespace xoos::cnc
