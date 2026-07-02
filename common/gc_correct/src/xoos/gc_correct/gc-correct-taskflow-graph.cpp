#include "xoos/gc_correct/gc-correct-taskflow-graph.h"

#include <taskflow/taskflow.hpp>

#include <xoos/log/logging.h>

#include "xoos/gc_correct/loess-taskflow-graph.h"

namespace xoos::gc_correct {

GCCorrectLoessInput GCCorrectTaskFlowGraph::PrepareDataForLoess(bool is_on_target) const {
  GCCorrectLoessInput res;
  res.is_valid_idxs = GetValidRegions(is_on_target);
  if (res.is_valid_idxs.empty()) {
    Logging::Info("no valid regions found for on-target={}. Skipping", is_on_target);
    return {};
  }
  res.is_ideal_idxs = GetIdealRegions(is_on_target, res.is_valid_idxs);
  res.ideal_average_coverage = _average_coverage.elem(res.is_ideal_idxs);
  res.ideal_gc_bias = _gc_bias.elem(res.is_ideal_idxs);
  res.valid_average_coverage = _average_coverage.elem(res.is_valid_idxs);
  res.valid_gc_bias = _gc_bias.elem(res.is_valid_idxs);
  return res;
}

void GCCorrectTaskFlowGraph::CorrectCoverageUsingSmoothedValues(const GCCorrectLoessInput& loess_input,
                                                                const arma::vec& smoothed_values) {
  arma::vec scale_factor = smoothed_values / arma::mean(loess_input.ideal_average_coverage);
  size_t j = 0;
  // Tracks whether the target is ideal or not.
  for (auto i : loess_input.is_valid_idxs) {
    if (scale_factor[j] > 0) {
      _res.total_coverage[i] = _total_coverage[i] / scale_factor[j];
      _res.counts[i] = _counts[i] / scale_factor[j];
      _res.average_coverage(i) = (1.0 * _res.total_coverage[i]) / static_cast<double>(_region_sizes[i]);
    } else {
      _res.total_coverage[i] = 0;
      _res.counts[i] = 0;
      _res.average_coverage(i) = 0;
    }
    j += 1;
  }
}

/**
 * @brief describe a 1-based closed genome interval. We see similar logic in
 * other parts of the codebase; I want to leave this as a private implementation
 * detail until I can figure out if there's a common struct that I can use
 * instead
 */
struct Region {
  std::string contig;
  arma::uword start;
  arma::uword end;
};

/**
 * @brief parses region string into contig, start and end. start is 0-based and end is 1-based exclusive (as in a BED
 * file)
 * @param r  region string in format "contig:start-end"
 * @return Region
 * Will also replace with a common function if one exists
 */
static Region ParseRegionString(const std::string& r) {
  auto colon_pos = r.find(':');
  auto dash_pos = r.find('-');
  auto contig = r.substr(0, colon_pos);
  auto start = std::stoul(r.substr(colon_pos + 1, dash_pos - colon_pos - 1)) - 1;
  auto end = std::stoul(r.substr(dash_pos + 1));
  return {.contig = contig, .start = start, .end = end};
}

arma::uvec GCCorrectTaskFlowGraph::GetRegionSizes(const std::vector<std::string>& regions) {
  vector<arma::uword> region_sizes;
  for (const std::string& region : regions) {
    Region reg = ParseRegionString(region);
    region_sizes.push_back(reg.end - reg.start);
  }
  return region_sizes;
}

/**
 * @brief using total coverages and region sizes, calculate the mean per-base coverage for each region
 * @return arma::vec of average coverages
 */
arma::vec GCCorrectTaskFlowGraph::CalculateAverageCoverage() const {
  vector<double> avg_cov;
  for (size_t i = 0; i < _region_sizes.size(); ++i) {
    arma::uword region_size = _region_sizes[i];
    avg_cov.push_back(static_cast<double>(_total_coverage[i]) / static_cast<double>(region_size));
  }
  return avg_cov;
}

/**
 * @brief contains the logic for performing GC correction, in the form of a taskflow graph
 * @return
 */
tf::Graph GCCorrectTaskFlowGraph::BuildGraph() {
  tf::Graph graph;
  tf::FlowBuilder builder(graph);
  tf::Task init_data_task = builder.emplace([this]() {
    _res.total_coverage = arma::vec(_total_coverage);
    _res.counts = arma::vec(_counts);
    _res.average_coverage = arma::vec(_average_coverage);
  });
  vector<tf::Task> gc_correct_tasks;
  gc_correct_tasks.reserve(2);
  // define the loess smoothing lambda
  for (bool is_on_target : {true, false}) {
    gc_correct_tasks.emplace_back(builder.emplace([this, is_on_target](tf::Subflow& subflow) {
      GCCorrectLoessInput loess_input{PrepareDataForLoess(is_on_target)};
      if (loess_input.is_valid_idxs.empty()) {
        return;
      }
      arma::vec gc_to_est = arma::regspace(kGCCorrectLoessGCMin, kGCCorrectLoessGCStep, kGCCorrectLoessGCMax);
      LoessTaskFlowGraph loess_round1_graph(
          gc_to_est, loess_input.ideal_gc_bias, loess_input.ideal_average_coverage, _first_span, kGCCorrectLoessDegree);
      LoessTaskFlowGraph loess_round2_graph(loess_input.valid_gc_bias,
                                            gc_to_est,
                                            loess_round1_graph.GetResult(),
                                            kGCCorrectSecondLoessSpan,
                                            kGCCorrectLoessDegree);
      tf::Task loess_task1 = subflow.composed_of(loess_round1_graph);
      tf::Task loess_task2 = subflow.composed_of(loess_round2_graph);
      tf::Task correction_task = subflow.emplace([&loess_input, &loess_round2_graph, this]() {
        CorrectCoverageUsingSmoothedValues(loess_input, loess_round2_graph.GetResult());
      });
      loess_task1.precede(loess_task2);
      loess_task2.precede(correction_task);
      subflow.join();
    }));
  }
  for (auto& gc_correct_task : gc_correct_tasks) {
    init_data_task.precede(gc_correct_task);
  }
  return graph;
}

size_t GCCorrectTaskFlowGraph::GetSize() const {
  return _counts.size();
}

/**
 * @brief gets "valid" regions - i.e. regions with valid coverage and GC bias values. Only these regions will be used
 * for GC correction
 * @param is_on_target - if true, restrict results to on-target regions. If false, restrict to off-target regions
 * @return
 */
arma::uvec GCCorrectTaskFlowGraph::GetValidRegions(bool is_on_target) const {
  vector<arma::uword> is_valid_idxs;
  for (size_t i = 0; i < _counts.size(); ++i) {
    if (_on_target[i] == is_on_target && (_average_coverage[i] > 0) && _gc_bias[i] >= 0) {
      is_valid_idxs.push_back(i);
    }
  }
  return {is_valid_idxs};
}

/**
 * @brief return valid (see definition of "valid" above) regions that are also "ideal" - i.e. regions where coverage and
 * GC bias are not outliers
 * @param is_on_target - if true, restrict results to on-target regions. If false, restrict to off-target regions
 * @param is_valid_idxs - list of indexes of "valid" regions
 * @return list of indexes of "ideal" regions
 */
arma::uvec GCCorrectTaskFlowGraph::GetIdealRegions(bool is_on_target, const arma::uvec& is_valid_idxs) const {
  vector<arma::uword> is_ideal_idxs;
  // calculate quantiles for avg coverage and gc bias
  arma::vec valid_average_coverage = _average_coverage.elem(is_valid_idxs);
  arma::vec valid_gc_bias = _gc_bias.elem(is_valid_idxs);
  // ranges for filtering out coverage and gc content outliers
  arma::vec likely_cov_range = arma::quantile(valid_average_coverage, arma::vec({0, 0.99}));
  arma::vec likely_gc_range = arma::quantile(valid_gc_bias, arma::vec({0.001, 1 - 0.001}));
  arma::uword likely_region_size = 0;
  // filter out outliers w.r.t region size for off-target regions
  if (!is_on_target) {
    arma::uvec valid_region_size = _region_sizes.elem(is_valid_idxs);
    arma::vec range = arma::quantile(valid_region_size, arma::vec({0.1}));
    likely_region_size = static_cast<arma::uword>(range[0]);
  }
  // gather "ideal" regions - i.e. valid regions w/o outliers
  for (auto i : is_valid_idxs) {
    if (((!is_on_target && _region_sizes[i] >= likely_region_size) || (is_on_target && _mappability[i] >= 1)) &&
        _is_in_autosome[i] && _average_coverage[i] >= likely_cov_range[0] &&
        _average_coverage[i] <= likely_cov_range[1] && _gc_bias[i] >= likely_gc_range[0] &&
        _gc_bias[i] <= likely_gc_range[1]) {
      is_ideal_idxs.push_back(i);
    }
  }
  return {is_ideal_idxs};
}

static const std::set<std::string> kAllosomes = {"chrY", "chrX"};

std::vector<bool> GCCorrectTaskFlowGraph::IsInAutosome(const std::vector<std::string>& regions) {
  std::vector<bool> ret(regions.size());
  for (size_t i = 0; i < regions.size(); ++i) {
    auto [contig, start, end] = ParseRegionString(regions[i]);
    ret[i] = (kAllosomes.find(contig) == kAllosomes.end());
  }
  return ret;
}

}  // namespace xoos::gc_correct
