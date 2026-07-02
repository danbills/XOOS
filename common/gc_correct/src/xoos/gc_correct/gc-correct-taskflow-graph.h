#pragma once

#include <armadillo>
#include <filesystem>
#include <string>
#include <vector>

#include <taskflow/core/graph.hpp>

#include "xoos/gc_correct/gc-correct-results.h"

using std::vector;
using Path = std::filesystem::path;

namespace xoos::gc_correct {

const double kGCCorrectSecondLoessSpan = 0.3;
const int kGCCorrectLoessDegree = 2;
const double kGCCorrectLoessGCMin = 0;
const double kGCCorrectLoessGCMax = 1;
const double kGCCorrectLoessGCStep = 0.001;

struct GCCorrectLoessInput {
  arma::uvec is_valid_idxs;
  arma::uvec is_ideal_idxs;
  arma::vec ideal_average_coverage;
  arma::vec ideal_gc_bias;
  arma::vec valid_average_coverage;
  arma::vec valid_gc_bias;
};

/**
 * @class  GCCorrectTaskFlowGraph
 * @brief performs GC bias correction on coverages
 */
class GCCorrectTaskFlowGraph {
 public:
  GCCorrectTaskFlowGraph() = delete;

  /**
   * @brief constructor -
   * @param counts
   * @param total_coverage
   * @param gc_bias
   * @param mappability
   * @param on_target
   * @param first_span - span for the first round of smoothing
   */
  GCCorrectTaskFlowGraph(const std::vector<std::string>& regions,
                         const arma::vec& counts,
                         const arma::vec& total_coverage,
                         const arma::vec& gc_bias,
                         const arma::vec& mappability,
                         const std::vector<bool>& on_target,
                         double first_span)
      : _counts(counts),
        _total_coverage(total_coverage),
        _gc_bias(gc_bias),
        _mappability(mappability),
        _on_target(on_target),
        _region_sizes(GetRegionSizes(regions)),
        _average_coverage(CalculateAverageCoverage()),
        _first_span(first_span),
        _is_in_autosome(IsInAutosome(regions)),
        _graph(BuildGraph()) {
  }

  const GCCorrectResults& GetResult() const {
    return _res;
  }

  GCCorrectResults& GetResult() {
    return _res;
  }

  /**
   * @brief gets the corrected AVERAGE coverages
   * @param on_target_only - specify whether to return only values in on-target
   *                         regions or not
   * @return  a list of the AVERAGE corrected coverages
   */
  size_t GetSize() const;

  // NOLINTNEXTLINE(readability-identifier-naming)
  tf::Graph& graph() {
    return _graph;
  }

 private:
  tf::Graph BuildGraph();
  arma::vec CalculateAverageCoverage() const;
  static arma::uvec GetRegionSizes(const std::vector<std::string>& regions);
  /**
   * @brief get indexes of all "valid" regions.
   * @param is_on_target - speficy whether to work on on-target region only
   *                       (true) or off-target regions only (false)
   * @return list of indexes for valid regions
   * A valid region is one where the coverage and GC content are valid numbers
   * (i.e. coverage >= 0 and GC content > 0)
   */
  arma::uvec GetValidRegions(bool is_on_target) const;
  /**
   * @brief get indexes of all "ideal" regions.
   * @param is_on_target - speficy whether to work on on-target region only
   *                       (true) or off-target regions only (false)
   * @param uvec - list of indexes for valid regions
   * An ideal region is one where the coverage falls within expected
   * boundaries. More specifically, it is a region that falls within the 0.1
   * and 0.9 quantiles of coverages over all regions specified
   * args: is_on_target: bool - if true, calculates ideal regions over on-target
   * regions only. If false, calculates ideal regions over off-target regions
   * only
   * @return list of indexes for ideal regions
   */
  arma::uvec GetIdealRegions(bool is_on_target, const arma::uvec&) const;
  GCCorrectLoessInput PrepareDataForLoess(bool is_on_target) const;
  void CorrectCoverageUsingSmoothedValues(const GCCorrectLoessInput& loess_input, const arma::vec& smoothed_values);
  static std::vector<bool> IsInAutosome(const std::vector<std::string>& regions);
  const arma::vec& _counts;
  const arma::vec& _total_coverage;
  const arma::vec& _gc_bias;
  const arma::vec& _mappability;
  const std::vector<bool>& _on_target;
  arma::uvec _region_sizes;
  arma::vec _average_coverage;
  GCCorrectResults _res;
  double _first_span;
  std::vector<bool> _is_in_autosome;
  tf::Graph _graph;
};
}  // namespace xoos::gc_correct
