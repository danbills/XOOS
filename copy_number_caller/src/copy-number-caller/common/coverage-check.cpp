#include "copy-number-caller/common/coverage-check.h"

namespace xoos::cnc {

CoverageCheckResult CheckMedianCoverage(const CoverageRecords& coverage, f64 threshold) {
  CoverageCheckResult result{};
  if (coverage.count.empty()) {
    result.median_coverage = 0.0;
    result.pct_windows_below_threshold = 100.0;
    result.is_sufficient = false;
    return result;
  }
  result.median_coverage = arma::median(coverage.count);
  const arma::uword n_below = arma::sum(coverage.count < threshold);
  result.pct_windows_below_threshold = 100.0 * static_cast<f64>(n_below) / static_cast<f64>(coverage.count.n_elem);
  result.is_sufficient = result.median_coverage >= threshold;
  return result;
}

}  // namespace xoos::cnc
