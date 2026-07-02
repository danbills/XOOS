#pragma once

#include <armadillo>

namespace xoos::gc_correct {
struct GCCorrectResults {
  arma::vec counts;
  arma::vec average_coverage;
  arma::vec total_coverage;
};
}  // namespace xoos::gc_correct
