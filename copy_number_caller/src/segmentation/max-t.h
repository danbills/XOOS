#pragma once

#include <armadillo>
#include <map>
#include <string>

#include "genomic-segments.h"

namespace xoos::cnc::segmentation {
enum class CbsMaxTMethod {
  kBruteForce,
  kFast,
  kFastMaxN
};
extern const std::map<std::string, CbsMaxTMethod> kStringToCbsMethod;
Segment GetMaxTBreakpoints(const arma::vec& obvs, size_t max_n, size_t min_n, CbsMaxTMethod method);
Segment GetMaxTBreakpointsBruteForce(const arma::vec& obvs, size_t max_n, size_t min_n);
Segment GetMaxTBreakpointsFast(const arma::vec& obvs, size_t max_n, size_t min_n);
Segment GetMaxTBreakpointsFastMaxN(const arma::vec& obvs, size_t max_n, size_t min_n);
}  // namespace xoos::cnc::segmentation
