#include "xoos/gc_correct/loess.h"

#include <taskflow/taskflow.hpp>

#include "xoos/gc_correct/loess-taskflow-graph.h"

namespace xoos::gc_correct {
arma::vec Loess(const arma::vec& x_to_predict,
                const arma::vec& x,
                const arma::vec& y,
                double first_span,
                int degree,
                int n_threads) {
  tf::Executor executor(n_threads);
  tf::Taskflow taskflow;
  LoessTaskFlowGraph loess_taskflow_graph(x_to_predict, x, y, first_span, degree);
  taskflow.composed_of(loess_taskflow_graph);
  executor.run(taskflow).get();
  return loess_taskflow_graph.GetResult();
}
}  // namespace xoos::gc_correct
