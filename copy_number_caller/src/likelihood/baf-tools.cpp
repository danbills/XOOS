#include "likelihood/baf-tools.h"

#include <cassert>

#include <xoos/log/logging.h>
#include <xoos/util/math.h>

namespace xoos::cnc {

/**
 * @brief calculates bandwidth for kernel density estimation using Silverman's rule of thumb.
  // https://en.wikipedia.org/wiki/Kernel_density_estimation#A_rule-of-thumb_bandwidth_estimator
 * @param x  data for which kernel density estimation is to be calculated
 * @return bandwidth parameter to use for kernel density estimation
 */
f64 GetBandwidthSilvermanRuleOfThumb(const arma::vec& x) {
  // Rule of thumb bandwidth estimator
  auto n = static_cast<f64>(x.n_elem);
  f64 std_dev = arma::stddev(x);
  arma::vec lower_upper_quartiles = arma::quantile(x, arma::vec({0.25, 0.75}));
  f64 iqr_coeff = (lower_upper_quartiles[1] - lower_upper_quartiles[0]) / 1.34;
  f64 a = std::min(std_dev, iqr_coeff);

  // Additional handling for a == 0
  if (math::IsCloseToZero(a)) {
    a = std_dev;
    if (math::IsCloseToZero(a)) {
      a = std::abs(x.n_elem > 0 ? x[0] : 0.0);
      if (math::IsCloseToZero(a)) {
        a = 1.0;
      }
    }
  }

  return 0.9 * a * std::pow(n, -1.0 / 5.0);
}

/**
 * @brief return the Guassian Kernel Density Estimate of x, weigted by `weights`. Bandwidth is calculated using the
 * rule-of-thumb bandwidth estimator
 * https://en.wikipedia.org/wiki/Kernel_density_estimation#A_rule-of-thumb_bandwidth_estimator
 * Binned KDE performed using procedure described in https://www.sfu.ca/sasdoc/sashtml/stat/chap33/sect12.htm
 * @param x data for which to calculate kernel density estimates
 * @param weights weights for each data in x
 * @param bandwidth_multiplier scalar applied to the Silverman rule-of-thumb bandwidth; default 1.0; must be > 0
 * @return KernelDensityEstimateRet -
 *          .x BAF bins at which density was estiamted. The value of each BAF bin is the weighted average of all the BAF
 * values in the same bin .y the density at each bin
 */
KernelDensityEstimateRet GaussianKernelDensityEstimate(const arma::vec& x, const arma::vec& weights) {
  return GaussianKernelDensityEstimate(x, weights, 1.0);
}

KernelDensityEstimateRet GaussianKernelDensityEstimate(const arma::vec& x,
                                                       const arma::vec& weights,
                                                       const f64 bandwidth_multiplier) {
  f64 weights_sum = arma::sum(weights);
  f64 diff_from_expected = fabs(weights_sum - static_cast<f64>(x.size()));
  if (!weights.empty() && diff_from_expected > 0.00001) {
    Logging::Error("GaussianKernelDensityEstimate: weights are expected to sum to x.size()!!!");
    throw std::runtime_error("weights do not sum to 1");
  }
  // Sort the x values and weights
  arma::uvec sorted_indices = arma::sort_index(x);
  arma::vec sorted_x = x.elem(sorted_indices);
  arma::vec sorted_weights = weights.empty()
                                 ? arma::vec(sorted_x.n_elem, arma::fill::ones) / static_cast<f64>(sorted_x.n_elem)
                                 : arma::vec(weights.elem(sorted_indices));
  // Bandwidth for Gaussian kernel
  const f64 bandwidth = GetBandwidthSilvermanRuleOfThumb(x) * bandwidth_multiplier;
  auto n = static_cast<f64>(x.n_elem);
  f64 coef = 1 / (n * bandwidth);

  // Calculate extension and new min/max for kde_x
  f64 extension = 3.0 * bandwidth;
  f64 min_x = sorted_x.min() - extension;
  f64 max_x = sorted_x.max() + extension;

  arma::vec density = arma::vec(512, arma::fill::zeros);
  arma::vec kde_x = arma::linspace(min_x, max_x, 512);
  for (size_t i = 0; i < kde_x.n_elem; ++i) {
    density(i) = arma::sum(sorted_weights % arma::normpdf((kde_x(i) - sorted_x) / bandwidth));
  }
  density = density * coef;
  return KernelDensityEstimateRet{.x = kde_x, .y = density};
}

/**
 * @brief calculate BAF values given ref and alt allelic depths
 * @param ref_ads  ref depths for a list of alleles
 * @param alt_ads  alt depths for the same list of alleles
 * @return std::vector of B-allele fractions (alt depth over total depth) - unmirrored.
 */
arma::vec GetBAFFromAD(const arma::vec& ref_ads, const arma::vec& alt_ads) {
  arma::vec bafs(ref_ads.size());
  // calculate BAF from ref_ads and alt_ads
  for (size_t i = 0; i < ref_ads.n_elem; ++i) {
    f64 total_depth = ref_ads(i) + alt_ads(i);
    if (math::IsCloseToZero(total_depth)) {
      bafs(i) = NAN;
    } else {
      bafs(i) = alt_ads(i) / total_depth;
    }
  }
  return bafs;
}

/**
 * @brief find indexes local maxima in a given array. If there's a "plateau" in the data where the left side of the
 * plateau is increasing and the right side is decreasing, count the median of the plateau as a local maxima.
 * @param x std::vector in which to find local maxima
 * @return std::vector of indices of local maxima
 */
std::vector<size_t> FindLocalMaxima(const arma::vec& x) {
  // state 0 = increasing
  // state 1 = no change (plateau)
  // state 2 = decreasing
  size_t state = 0;
  std::vector<size_t> local_maxima_idxs;
  if (x[0] > x[1]) {
    local_maxima_idxs.push_back(0);
    state = 2;
  }
  size_t plateau_start_idx = 0;
  bool plateau_possibly_maximum = false;
  for (size_t i = 1; i < x.size() - 1; ++i) {
    if (state == 0) {  // previously increasing
      if (math::IsEqualWithTolerance(x[i], x[i - 1])) {
        state = 1;
        plateau_start_idx = i - 1;
        plateau_possibly_maximum = true;
      } else if (x[i] < x[i - 1]) {  // immediately known to be a local maximum
        state = 2;
        local_maxima_idxs.push_back(i - 1);
      }  // otherwise keep state at 0 and do nothing
    } else if (state == 1) {  // previously "equal"
      if (x[i] < x[i - 1]) {
        state = 2;  // now decreasing, we can set the middle of the plateau as a local maximum
        if (plateau_possibly_maximum) {
          size_t idx = plateau_start_idx + (i - plateau_start_idx) / 2;
          local_maxima_idxs.push_back(idx);
        }
      } else if (x[i] > x[i - 1]) {
        state = 0;  // do nothing else since this is not a maximum
      }
    } else if (state == 2) {
      if (x[i] > x[i - 1]) {
        state = 0;  // now increasing, but we don't care about minima here so do nothing
      } else if (math::IsEqualWithTolerance(x[i], x[i - 1])) {
        state = 1;  // plateau, do nothing
        plateau_possibly_maximum = false;
      }  // otherwise keep state at 2 and do nothing
    }
  }
  return local_maxima_idxs;
}

/**
 * @brief find the peak of a BAF distribution. This is meant as an alternative
 * to getting the mean BAF after mirroring. By finding the peak of the
 * distribution and then mirroring, we provide a more accurate estimation of the
 * mean mirrored baf
 * @param ref_ads
 * @param alt_ads
 * @return
 */
f64 GetPeakOfBAFDistributionThenMirror(const arma::vec& ref_ads, const arma::vec& alt_ads) {
  arma::vec bafs = GetBAFFromAD(ref_ads, alt_ads);
  // add in 1-baf in order to f64 the data for calculating density (pronounce the curves more)
  arma::vec aafs = 1 - bafs;
  arma::vec kde_input = arma::join_cols(bafs, aafs);
  arma::uvec sort_order = arma::sort_index(kde_input);
  // weight the aaf+baf observations by how many reads overlap that observation
  arma::vec total_depth = ref_ads + alt_ads;
  total_depth = arma::join_cols(total_depth, total_depth);
  arma::vec weights = total_depth / (arma::sum(total_depth)) * static_cast<f64>(total_depth.n_elem);
  // sort all the input
  kde_input = kde_input.elem(sort_order);
  weights = weights.elem(sort_order);
  total_depth = total_depth.elem(sort_order);
  auto kde = GaussianKernelDensityEstimate(kde_input, weights);
  // for now, just simply calculate all local maxima, and then take the max of those.
  // if needed, we could find "strong" maxima, but at this moment I think that is overkill
  auto local_maxima = FindLocalMaxima(kde.y);
  // If failed to find local maxima
  if (local_maxima.empty()) {
    local_maxima.push_back(kde.y.index_max());
  }
  // we should be suspicious of any segment that has more than 2 peaks - could indicate a segment that could be broken
  // down even further
  if (local_maxima.size() > 2) {
    Logging::Debug("Segment has more than 2 peaks in its BAF distribution");
  }
  size_t global_maximum_idx = 0;
  for (auto local_maximum_idx : local_maxima) {
    if (kde.y[local_maximum_idx] > kde.y[global_maximum_idx]) {
      global_maximum_idx = local_maximum_idx;
    }
  }
  f64 max_mirrored_baf = kde.x[global_maximum_idx];
  if (max_mirrored_baf > 0.5) {
    // if the maximum is above 0.5, we should mirror it
    max_mirrored_baf = 1 - max_mirrored_baf;
  }
  return max_mirrored_baf;
}

f64 WeightedMBAFSD(const arma::vec& ref_ads, const arma::vec& alt_ads, f64 mean_baf) {
  f64 num = 0;
  f64 denom = 0;
  for (size_t i = 0; i < ref_ads.size(); ++i) {
    f64 minor_depth = std::min({ref_ads[i], alt_ads[i]});
    f64 total_depth = ref_ads[i] + alt_ads[i];
    if (minor_depth < 0) {
      Logging::Error("Found a reference or alternative allele depth that is < 0");
      throw std::runtime_error("Found a reference or alternative allele depth that is < 0");
    }
    num += sqrt(total_depth) * pow((minor_depth / total_depth) - mean_baf, 2);
    denom += sqrt(total_depth);
  }
  return sqrt(num / (denom - 1));
}

}  // namespace xoos::cnc
