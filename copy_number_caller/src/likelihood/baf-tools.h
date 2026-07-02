#pragma once
#include <armadillo>

#include <xoos/types/float.h>

namespace xoos::cnc {

struct KernelDensityEstimateRet {
  arma::vec x;
  arma::vec y;
};

KernelDensityEstimateRet GaussianKernelDensityEstimate(const arma::vec& x, const arma::vec& weights);
KernelDensityEstimateRet GaussianKernelDensityEstimate(const arma::vec& x,
                                                       const arma::vec& weights,
                                                       f64 bandwidth_multiplier);
arma::vec GetBAFFromAD(const arma::vec& ref_ads, const arma::vec& alt_ads);
std::vector<size_t> FindLocalMaxima(const arma::vec& x);
f64 GetPeakOfBAFDistributionThenMirror(const arma::vec& ref_ads, const arma::vec& alt_ads);
f64 WeightedMBAFSD(const arma::vec& ref_ads, const arma::vec& alt_ads, f64 mean_baf);
f64 GetBandwidthSilvermanRuleOfThumb(const arma::vec& x);
}  // namespace xoos::cnc
