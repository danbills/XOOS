#pragma once
#include <string>
#include <unordered_map>
#include <vector>

namespace xoos::cnc {

struct Interval {
  size_t start, end;
  auto operator<=>(const Interval&) const = default;
};

using ContigToIntervals = std::unordered_map<std::string, std::vector<Interval>>;

};  // namespace xoos::cnc
