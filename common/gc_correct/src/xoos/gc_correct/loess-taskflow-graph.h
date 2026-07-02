#pragma once
#include <armadillo>

#include <taskflow/core/graph.hpp>

namespace xoos::gc_correct {
class LoessTaskFlowGraph {
 public:
  LoessTaskFlowGraph() = delete;

  LoessTaskFlowGraph(const arma::vec& x_to_predict, const arma::vec& x, const arma::vec& y, double span, int degree)
      : _x_to_predict(x_to_predict), _x(x), _y(y), _span(span), _degree(degree), _graph(BuildGraph()) {
  }

  const arma::vec& GetResult() {
    return _ret;
  }

  // NOLINTNEXTLINE(readability-identifier-naming)
  tf::Graph& graph() {
    return _graph;
  }

 private:
  tf::Graph BuildGraph();
  double SmoothPoint(int i) const;
  const arma::vec& _x_to_predict;
  const arma::vec& _x;
  const arma::vec& _y;
  double _span;
  int _degree;
  arma::uvec _sort_order;
  arma::vec _sorted_xx;
  arma::vec _sorted_yy;
  int _n = 0;
  arma::vec _ret;
  tf::Graph _graph;
};
}  // namespace xoos::gc_correct
