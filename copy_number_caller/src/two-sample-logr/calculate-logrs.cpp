#include "two-sample-logr/calculate-logrs.h"

#include <armadillo>

#include <xoos/error/error.h>
#include <xoos/util/math.h>

#include "utility/utility-functions.h"

namespace xoos::cnc {
void FilterAndNormalizeTumorNormal(CoverageRecords& tumor,
                                   CoverageRecords& normal,
                                   size_t normal_min_cov,
                                   size_t tumor_min_cov) {
  std::vector<arma::uword> keep_idxs;
  for (size_t i = 0; i < normal.average_coverage.size(); ++i) {
    if (tumor.average_coverage[i] >= static_cast<f64>(tumor_min_cov) &&
        normal.average_coverage[i] >= static_cast<f64>(normal_min_cov)) {
      keep_idxs.push_back(i);
    }
  }
  if (keep_idxs.empty()) {
    throw error::Error("No regions passed coverage thresholds (normal >= {}, tumor >= {}); cannot normalize",
                       normal_min_cov,
                       tumor_min_cov);
  }
  normal = normal.FilterRow(keep_idxs);
  tumor = tumor.FilterRow(keep_idxs);
  const auto normal_sum = arma::sum(normal.count);
  const auto tumor_sum = arma::sum(tumor.count);
  if (math::IsCloseToZero(normal_sum) || math::IsCloseToZero(tumor_sum)) {
    throw error::Error(
        "Total count is zero after filtering (normal sum: {}, tumor sum: {}); cannot normalize", normal_sum, tumor_sum);
  }
  normal.count = normal.count / normal_sum;
  tumor.count = tumor.count / tumor_sum;
}

Observations CalculateLogrs(const CoverageRecords& tumor, const CoverageRecords& normal_normalized) {
  CoverageRecords tumor_normalized = tumor.FilterRegion(normal_normalized.region);
  Observations res;
  res.regions.resize(normal_normalized.region.size());
  res.contigs.resize(normal_normalized.region.size());
  res.starts.resize(normal_normalized.region.size());
  res.ends.resize(normal_normalized.region.size());
  res.obvs.resize(normal_normalized.region.size());
  for (size_t i = 0; i < normal_normalized.region.size(); ++i) {
    if (normal_normalized.region[i] != tumor_normalized.region[i]) {
      throw error::Error("Tumor and normal regions do not match at index {}: '{}' vs '{}'",
                         i,
                         tumor_normalized.region[i],
                         normal_normalized.region[i]);
    }
    const auto& [contig, start, end] = ParseRegionString(normal_normalized.region[i]);
    res.contigs[i] = contig;
    res.starts[i] = start;
    res.ends[i] = end;
    res.regions[i] = normal_normalized.region[i];
  }
  // calculate log2-ratio
  res.obvs = arma::log2(tumor_normalized.count) - arma::log2(normal_normalized.count);

  // replace log2-ratio values that are -inf with the minimum log2-ratio value
  // -inf can happen if the tumour coverage is 0X
  arma::uvec finite_indices = arma::find_finite(res.obvs);
  arma::vec finite_values(finite_indices.n_elem);
  finite_values = res.obvs.elem(finite_indices);

  for (auto& val : res.obvs) {
    if (std::isinf(val) && val < 0) {
      val = arma::min(finite_values);
    }
  }
  // median-center the log2-ratios
  res.obvs = res.obvs - arma::median(res.obvs);
  return res;
}
}  // namespace xoos::cnc
