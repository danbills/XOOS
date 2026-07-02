#include "filter-observations.h"

#include <vector>

#include "observations.h"
#include "segmentation/interval-trees.h"

namespace xoos::cnc {
using segmentation::IntervalTrees;

/**
 * @brief Get observations from `right` that overlap with `left`
 * @param left Observations object to filter from
 * @param right Observations object to filter with
 * @return Observations object containing observations from right that overlap with left. Sets IsBAFSegObvs
 * automatically
 */
Observations GetOverlappingObservations(const Observations& left, const Observations& right) {
  if (right.contigs.empty()) {
    return {};
  }
  Observations ret;
  std::vector<arma::uword> starts;
  std::vector<arma::uword> ends;
  std::vector<f64> obvs;
  std::vector<s32> dps;
  IntervalTrees left_trees(left);
  for (size_t i = 0; i < right.contigs.size(); ++i) {
    std::vector<size_t> idxs = left_trees.LookUp(right.contigs[i], right.starts[i], right.ends[i]);
    // we don't want to keep the observation if there are no overlapping regions
    if (idxs.empty()) {
      continue;
    } else {
      ret.contigs.push_back(right.contigs[i]);
      starts.push_back(right.starts[i]);
      ends.push_back(right.ends[i]);
      obvs.push_back(right.obvs[i]);
      if (!right.dps.empty()) {
        ret.dps.push_back(right.dps[i]);
      }
    }
  }
  ret.starts = arma::uvec(starts);
  ret.ends = arma::uvec(ends);
  ret.obvs = arma::vec(obvs);
  return ret;
}

/**
 * @brief Get observations from `right` that overlap with `left`
 * @param left list of region strings to filter from
 * @param right Observations object to filter with
 * @return Observations object containing observations from right that overlap with left. Does not set IsBAFSegObvs
 * automatically
 */
Observations GetOverlappingObservations(const std::vector<std::string>& left, const Observations& right) {
  if (right.contigs.empty()) {
    return {};
  }
  Observations ret;
  std::vector<arma::uword> starts;
  std::vector<arma::uword> ends;
  std::vector<f64> obvs;
  std::vector<s32> dps;
  IntervalTrees left_trees(left);
  for (size_t i = 0; i < right.contigs.size(); ++i) {
    std::vector<size_t> idxs = left_trees.LookUp(right.contigs[i], right.starts[i], right.ends[i]);
    // we don't want to keep the observation if there are no overlapping regions
    if (idxs.empty()) {
      continue;
    } else {  // keep the observation from `right`
      ret.contigs.push_back(right.contigs[i]);
      starts.push_back(right.starts[i]);
      ends.push_back(right.ends[i]);
      obvs.push_back(right.obvs[i]);
      if (!right.dps.empty()) {
        ret.dps.push_back(right.dps[i]);
      }
    }
  }
  ret.starts = arma::uvec(starts);
  ret.ends = arma::uvec(ends);
  ret.obvs = arma::vec(obvs);
  return ret;
}
}  // namespace xoos::cnc
