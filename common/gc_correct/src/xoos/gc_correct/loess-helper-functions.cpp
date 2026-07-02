#include "xoos/gc_correct/loess-helper-functions.h"

#include <xoos/log/logging.h>

namespace xoos::gc_correct {

/**
 * @brief Get the range of the top N closest elements to x in arr. Assume arr is sorted
 * @param arr
 * @param n
 * @param x
 * @return the first and last indices in arr representing N closest elements in arr to x. [first, last] is inclusive,
 * closed
 */
std::tuple<size_t, size_t> GetNClosestRange(const arma::vec& arr, size_t n, double x) {
  // since arr is sorted, arr - x will be sorted too, but abs(arr - x) will not be sorted. Will have to find location of
  // min, then do a search from there
  arma::uword first = 0;
  arma::uword last = arr.n_elem - 1;
  while (last - first >= n) {
    // if the distance`first` is greater than the distance at `last`, then we can move `first` up
    if (std::abs(arr[first] - x) > std::abs(arr[last] - x)) {
      first++;
    } else {
      last--;
    }
  }
  return {first, last};
}

arma::vec Tricubic(const arma::vec& arr) {
  return arma::pow(1 - arma::pow(arr, 3), 3);
}

arma::vec TricubicDistWeights(const arma::vec& v, double x) {
  const arma::vec distances = arma::abs(v - x);
  const double max_distance = arma::max(distances);
  // if max distance is 0, treat as if all normed distances are 1 away.
  if (max_distance == 0) {
    return {v.n_elem, arma::fill::ones};
  }
  const arma::vec normed_distances = distances / max_distance;
  return Tricubic(normed_distances);
}

arma::vec WeightedPolyFit(
    const arma::vec& x, const arma::vec& y, const arma::vec& new_x, const arma::vec& weights, arma::uword degree) {
  const arma::mat x_poly = MakePoly(x, degree);
  arma::vec beta(degree + 1);
  if (weights.n_elem != 0) {
    // LS estimate is: INV(T(X)WX)T(X)WY
    // Eq 7.4 from https://www.stat.cmu.edu/~cshalizi/uADA/16/lectures/08.pdf
    arma::mat xt_w = x_poly.t();
    // multiply x_poly_t by diag(w) (diagonal of weights) 3xn * nxn -> 3xn
    for (size_t i = 0; i < weights.n_elem; ++i) {
      xt_w.col(i) *= weights[i];
    }
    // solve: xt_w_x * beta = xt_w_y
    const arma::mat xt_w_x = xt_w * x_poly;
    const arma::mat xt_w_y = xt_w * y;
    beta = arma::solve(xt_w_x, xt_w_y, arma::solve_opts::fast);
  } else {
    // 3 * n
    const arma::mat xt = x_poly.t();
    // 3x3
    const arma::mat xt_x_inv = arma::pinv(xt * x_poly);
    // 3x3 * 3xn -> 3xn
    const arma::mat xt_x_inv_xt = xt_x_inv * xt;
    // 3xn * nx1 -> 3x1
    beta = xt_x_inv_xt * y;
  }
  // nx3
  const arma::mat poly_x = MakePoly(new_x, degree);
  // nx3 * 3x1 -> nx1
  return poly_x * beta;
}

arma::mat MakePoly(const arma::vec& x, size_t degree) {
  arma::mat x_poly(x.n_elem, degree + 1, arma::fill::none);
  // fill each column (idx=c) with successive polynomials d of x, 0<=d<=2
  for (size_t i = 0; i < degree + 1; ++i) {
    arma::vec to_insert = arma::pow(x, static_cast<double>(i));
    for (size_t j = 0; j < to_insert.n_elem; ++j) {
      x_poly(j, i) = to_insert[j];
    }
  }
  return x_poly;
}
}  // namespace xoos::gc_correct
