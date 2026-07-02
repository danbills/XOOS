#include "histogram-summary.h"

namespace xoos::histogram {

std::optional<f64> GetRatioPercentile(const PercentileMap& percentiles, const u64 numerator, const u64 denominator) {
  // check if the numerator and denominator are valid percentiles
  if (!percentiles.contains(numerator) || !percentiles.contains(denominator) ||
      !percentiles.at(numerator).has_value() || !percentiles.at(denominator).has_value()) {
    return std::nullopt;
  }
  // if the denominator is zero, return std::nullopt to avoid division by zero
  if (percentiles.at(denominator).value() == 0) {
    return std::nullopt;
  }
  return static_cast<f64>(percentiles.at(numerator).value()) / static_cast<f64>(percentiles.at(denominator).value());
}

bool OmitsFirstBin(const HistogramBinOutput mode) {
  return mode == kOmitFirstBin || mode == kOmitFirstBinAndMaxLastBinWithOutlier || mode == kOmitFirstBinAndMaxLastBin;
}

bool FoldsOutliersIntoLastBin(const HistogramBinOutput mode) {
  return mode == kMaxLastBin || mode == kOmitFirstBinAndMaxLastBin;
}

}  // namespace xoos::histogram
