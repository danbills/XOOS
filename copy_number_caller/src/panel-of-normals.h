#pragma once

#include <armadillo>
#include <string>
#include <vector>

#include <xoos/types/float.h>
#include <xoos/types/fs.h>

#include "coverage.h"
#include "observations.h"
#include "singular-value-decomposition.h"

namespace xoos::cnc {
const std::string kOnTargetPrefix = "on_target";
const std::string kOffTargetPrefix = "off_target";

class PanelOfNormals {
 public:
  PanelOfNormals(const std::vector<CoverageRecords>& covs, bool on_target);
  PanelOfNormals(const fs::path& fname, bool on_target);
  PanelOfNormals() = default;
  PanelOfNormals(const PanelOfNormals& rhs) = default;
  PanelOfNormals(PanelOfNormals&& rhs) noexcept = default;
  PanelOfNormals& operator=(const PanelOfNormals& rhs) = default;
  PanelOfNormals& operator=(PanelOfNormals&& rhs) noexcept = default;
  ~PanelOfNormals() noexcept = default;
  void WriteStandardizedCounts(std::ostream& ofs) const;
  void WriteOriginalCounts(std::ostream& ofs) const;
  void WriteNormedIntervalMedians(std::ostream& ofs) const;
  void WriteOriginalIntervalMedians(std::ostream& ofs) const;
  void WriteSex(std::ostream&) const;
  std::vector<std::string> GetIntervals() const;
  arma::vec GetNormedMaleIntervalMedians() const;
  arma::vec GetNormedFemaleIntervalMedians() const;
  arma::vec GetNormedAllIntervalMedians() const;
  arma::vec GetNormedSexIntervalMedians(Sex sex) const;
  arma::vec GetOriginalSexIntervalMedians(Sex sex) const;
  const SingularValueDecomposition& GetSingularValueDecomposition() const;
  Observations LogRatio(const CoverageRecords& sample, Sex sex) const;
  Observations DenoiseLogR(const Observations& log_ratios, arma::uword num_eigen) const;
  std::vector<std::string> FilterIntervalsFromSample(const CoverageRecords& sample,
                                                     Sex sex,
                                                     size_t min_target_length,
                                                     f64 min_panel_median_cov,
                                                     f64 min_panel_median_and_tumor_cov) const;

  const arma::mat& GetOriginalCounts() const {
    return _original_counts;
  }

  const arma::vec& GetOriginalIntervalMedians() const {
    return _original_interval_medians;
  }

  // NOLINTNEXTLINE
  bool empty() const {
    return _counts.empty();
  }

 private:
  void LoadReferenceCovs(const std::vector<CoverageRecords>& covs, bool on_target);
  void FilterCounts();
  static arma::uvec IdentifyIntervalsNonLowCov(const arma::mat& counts);
  static arma::uvec IdentifySamplesNonExtremeCov(const arma::mat& counts);
  static arma::uvec IdentifyIntervalsNonZeroCov(const arma::mat& counts);
  static arma::uvec IdentifyIntervalsNonZeroMedian(const arma::mat& counts);
  void SetIntervalMedians(const arma::mat& counts,
                          arma::vec& all_medians,
                          arma::vec& male_medians,
                          arma::vec& female_medians);
  void StandardizeReferenceCovs();
  static void ImputeCounts(arma::mat& normed_counts);
  static void TruncateOutlierCounts(arma::mat& normed_counts, f64 lo, f64 hi);
  void FilterIntvs(const arma::uvec& idxs);
  void FilterSamples(const arma::uvec& idxs);
  void CalculateSingularValueDecomposition();
  bool _on_target = true;
  arma::mat _original_counts;  // un-normalized counts
  arma::mat _counts;           // normalized, standardized, **log2-transformed** and median-centered counts
  arma::vec _original_interval_medians;
  arma::vec _original_male_interval_medians;
  arma::vec _original_female_interval_medians;
  arma::vec _normed_interval_medians;
  arma::vec _normed_male_interval_medians;
  arma::vec _normed_female_interval_medians;
  std::vector<std::string> _kept_names;
  std::vector<std::string> _kept_intvs;
  std::vector<Sex> _sex;
  std::vector<f64> _xy_ratios;
  SingularValueDecomposition _svd;
};
}  // namespace xoos::cnc
