#pragma once
#include <armadillo>
#include <string>
#include <vector>

#include "xoos/gc_correct/gc-correct-results.h"

namespace xoos::gc_correct {

GCCorrectResults GCCorrect(const std::vector<std::string>& regions,
                           const arma::vec& counts,
                           const arma::vec& total_coverage,
                           const arma::vec& gc_bias,
                           const arma::vec& mappability,
                           const std::vector<bool>& on_target,
                           int n_threads,
                           double first_span);
}  // namespace xoos::gc_correct
