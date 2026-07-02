#pragma once
#include <armadillo>

#include <xoos/types/float.h>

namespace xoos::cnc::segmentation {
f64 TStat(f64 n_v1, f64 mean_v1, f64 var_v1, f64 n_v2, f64 mean_v2, f64 var_v2);
f64 TStatFromVecs(const arma::vec& v1, const arma::vec& v2);
f64 Variance(f64 m, f64 mean, f64 sq_sum);
arma::vec& FisherYatesShuffleInPlace(arma::vec& arr);
f64 TStatSquaredMean(const arma::vec& psums, size_t i, size_t j);
f64 TStatSquaredMean(f64 total_sum, f64 total_length, f64 minor_sum, f64 minor_length);
f64 TStatSquaredMean(f64 total_sum, size_t total_length, f64 minor_sum, size_t minor_length);
}  // namespace xoos::cnc::segmentation
