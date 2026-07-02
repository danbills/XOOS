#pragma once

#include <armadillo>

namespace xoos::gc_correct {
arma::vec Loess(const arma::vec& x_to_predict,
                const arma::vec& x,
                const arma::vec& y,
                double first_span,
                int degree,
                int n_threads);
}
