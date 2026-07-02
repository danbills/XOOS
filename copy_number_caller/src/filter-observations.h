#pragma once
#include <string>
#include <vector>

#include "observations.h"

namespace xoos::cnc {
Observations GetOverlappingObservations(const Observations& left, const Observations& right);
Observations GetOverlappingObservations(const std::vector<std::string>& left, const Observations& right);
}  // namespace xoos::cnc
