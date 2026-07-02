#include "xoos/gc_correct/gc-correct.h"

#include <taskflow/taskflow.hpp>

#include "xoos/gc_correct/gc-correct-taskflow-graph.h"

namespace xoos::gc_correct {
/**
 * @brief perform gc correction, given coverages and appropriate metadata
 * @param regions region strings
 * @param counts  read counts for each region
 * @param total_coverage total coverage (bases) for each region
 * @param gc_bias  gc bias fraction for each region
 * @param mappability  mappability score for each region
 * @param on_target <true> if "on" target, and <false> if "off" target. Set to <true> for all regions if not known.
 * @param n_threads number of threads
 * @param first_span suggested value: 0.03. To be deprecated
 * @return GCCorrectResults object (see gc-correct-results.h)
 */
GCCorrectResults GCCorrect(const std::vector<std::string>& regions,
                           const arma::vec& counts,
                           const arma::vec& total_coverage,
                           const arma::vec& gc_bias,
                           const arma::vec& mappability,
                           const std::vector<bool>& on_target,
                           int n_threads,
                           double first_span) {
  tf::Taskflow taskflow;
  tf::Executor executor(n_threads);
  GCCorrectTaskFlowGraph gc_correction_taskflow_graph(
      regions, counts, total_coverage, gc_bias, mappability, on_target, first_span);
  taskflow.composed_of(gc_correction_taskflow_graph);  // NOLINT
  executor.run(taskflow).get();
  GCCorrectResults res = gc_correction_taskflow_graph.GetResult();
  return res;
}
}  // namespace xoos::gc_correct
