#pragma once

#include <armadillo>
#include <tuple>

namespace xoos::gc_correct {

/**
 *  @brief return range of n closest indexes in x to index i, including i, where
 *        distance = abs(x[j] - x[i])
 *
 * arguments: sorted list, desired number of closest neighbors to find, index of
 * value whose neighbors need to be found
 *
 */
std::tuple<size_t, size_t> GetNClosestRange(const arma::vec& arr, size_t n, double x);

/**
 * @brief performs tricubic weight function on distances of v to x
 */
arma::vec TricubicDistWeights(const arma::vec& v, double x);

/** @brief tricubic weight function (1-d^3)^3 */
arma::vec Tricubic(const arma::vec&);

/**
 * @brief weighted polynomial fit
 * @param x  - predictor vector
 * @param y  - response vector
 * @param new_x - x to predict on
 * @param weights  - arbitrary weights per x value
 * @param degree  - polynomial degree
 * @return  new response vector representing values on fitted curve
 */
arma::vec WeightedPolyFit(
    const arma::vec& x, const arma::vec& y, const arma::vec& new_x, const arma::vec& weights, arma::uword degree);

/**
 * @brief Make a Vendermonde matrix for v
 */
arma::mat MakePoly(const arma::vec& x, size_t degree);

}  // namespace xoos::gc_correct
