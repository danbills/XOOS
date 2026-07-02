#include "segmentation/cbs-util.h"

#include <cmath>

#include <xoos/types/int.h>

namespace xoos::cnc::segmentation {
/**
 * @brief Imagine segment as circle joined at both ends. This function calculates t-statistic for testing the hypothesis
 * that arc from i+1 to j and its complement have different means.
 * @param n_v1 size of minor arc
 * @param mean_v1 mean of minor arc
 * @param var_v1 variance of minor arc
 * @param n_v2 size of major arc
 * @param mean_v2 mean of major arc
 * @param var_v2 variance of major arc
 * @return t-statistic describing difference in means between minor and major arc
 */
f64 TStat(f64 n_v1, f64 mean_v1, f64 var_v1, f64 n_v2, f64 mean_v2, f64 var_v2) {
  // pooled variance of v1 and v2
  f64 s_ij = std::pow(((n_v1 - 1) * var_v1 + (n_v2 - 1) * var_v2) / (n_v1 + n_v2 - 2), 0.5);
  // t stat is a function of the means of both std::vectors and their pooled variance
  f64 t_ij = (mean_v1 - mean_v2) / (s_ij * std::pow((1.0 / n_v1 + 1.0 / n_v2), 0.5));
  return t_ij;
}

/**
 * @brief calculates t-statistic between two std::vectors that can be of unequal size and variance
 * @param v1 std::vector 1
 * @param v2 std::vector 2
 * @return t-statistic between v1 and v2
 */
f64 TStatFromVecs(const arma::vec& v1, const arma::vec& v2) {
  auto n_v1 = static_cast<f64>(v1.n_elem);
  auto n_v2 = static_cast<f64>(v2.n_elem);
  f64 mean_v1 = arma::mean(v1);
  f64 mean_v2 = arma::mean(v2);
  f64 var_v1 = arma::var(v1);
  f64 var_v2 = arma::var(v2);
  return TStat(n_v1, mean_v1, var_v1, n_v2, mean_v2, var_v2);
}

/* warning - this is the biased estimator for variance. Do we need an unbiased estimator??? */
f64 Variance(f64 m, f64 mean, f64 sq_sum) {
  return (sq_sum / m) - std::pow(mean, 2);
}

arma::vec& FisherYatesShuffleInPlace(arma::vec& arr) {
  for (s32 i = static_cast<s32>(arr.size()) - 1; i > 0; --i) {
    // get a random number
    size_t random_index = arma::randi(arma::distr_param(0, i));
    std::swap(arr[i], arr[random_index]);
  }
  return arr;
}

/**
 * @brief return a t-stat that the difference in means between a "within" sub-segment (specified by i,j) against its
 * complementary "without" segment
 *
 * sub-segmnets. this generally works if we're looking for small sub-segments within very large segments
 * @param psums
 * @param i
 * @param j
 * @return
 */
f64 TStatSquaredMean(const arma::vec& psums, size_t i, size_t j) {
  auto inner_length = static_cast<f64>(j - i);
  auto outer_length = static_cast<f64>(psums.size() - j + i);
  f64 sum1 = i > 0 ? psums[j - 1] - psums[i - 1] : psums[j - 1];
  f64 sum2 = psums[psums.size() - 1] - sum1;
  return std::pow(sum1 / inner_length - sum2 / outer_length, 2) / (1 / inner_length + 1 / outer_length);
}

f64 TStatSquaredMean(f64 total_sum, f64 total_length, f64 minor_sum, f64 minor_length) {
  f64 major_sum = total_sum - minor_sum;
  f64 major_length = total_length - minor_length;
  return std::pow(major_sum / major_length - minor_sum / minor_length, 2) / (1 / major_length + 1 / minor_length);
}

f64 TStatSquaredMean(f64 total_sum, size_t total_length, f64 minor_sum, size_t minor_length) {
  return TStatSquaredMean(total_sum, static_cast<f64>(total_length), minor_sum, static_cast<f64>(minor_length));
}

}  // namespace xoos::cnc::segmentation
