#include "xoos/gc_correct/loess-taskflow-graph.h"

#include <taskflow/algorithm/for_each.hpp>

#include "xoos/gc_correct/loess-helper-functions.h"

namespace xoos::gc_correct {
double LoessTaskFlowGraph::SmoothPoint(int i) const {
  // does sorting help here?
  auto [start, end] = GetNClosestRange(_sorted_xx, _n, _x_to_predict[i]);
  if (_n > 1) {
    const arma::vec xx_subvec = _sorted_xx.subvec(arma::span(start, end));
    const arma::vec yy_subvec = _sorted_yy.subvec(arma::span(start, end));
    const arma::vec weights = TricubicDistWeights(xx_subvec, _x_to_predict[i]);
    const double y_pred = WeightedPolyFit(xx_subvec, yy_subvec, arma::vec({_x_to_predict[i]}), weights, _degree)[0];
    return y_pred;
  } else {
    return _sorted_yy[start];
  }
}

tf::Graph LoessTaskFlowGraph::BuildGraph() {
  tf::Graph graph;
  tf::FlowBuilder builder(graph);
  tf::Task setup = builder.emplace([this]() {
    // prepare the input data
    _sort_order = arma::sort_index(_x);
    _sorted_xx = _x.elem(_sort_order);
    _sorted_yy = _y.elem(_sort_order);
    _n = static_cast<int>(_span * static_cast<double>(_x.n_elem));
    if (_n == 0) {
      throw std::runtime_error("Too few data points per window. Try increasing the --first-span parameter.");
    }
    _ret.resize(_x_to_predict.n_elem);
  });
  const tf::Task predict =
      builder.for_each_index(0, static_cast<int>(_x_to_predict.n_elem), 1, [this](int i) { _ret[i] = SmoothPoint(i); });
  setup.precede(predict);
  return graph;
}
}  // namespace xoos::gc_correct
