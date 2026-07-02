#include "panel-of-normals.h"

#include <cassert>
#include <vector>

#include <xoos/log/logging.h>
#include <xoos/util/math.h>

#include "sex.h"
#include "utility/utility-functions.h"

namespace xoos::cnc {

const f64 kSampleOutlierMinPct = 0.25;
const f64 kSampleOutlierMaxPct = 4;
const f64 kIntvOutlierMinZeroPct = 0.03;
const f64 kIntvOutlierMin = 0.25;

/**
 * @brief load a panel file containing paths of reference coverage files,
 * standardize, and perform SVD
 * @param pon_fname - path name of reference panel
 * @param on_target - load on-target intervals only if true, else off-targets only
 */
PanelOfNormals::PanelOfNormals(const std::vector<CoverageRecords>& covs, bool on_target) {
  _on_target = on_target;
  Logging::Info("loading the reference coverages");
  LoadReferenceCovs(covs, on_target);
  if (!_counts.empty()) {
    Logging::Info("filtering counts");
    FilterCounts();
    Logging::Info("standardizing counts");
    StandardizeReferenceCovs();
    Logging::Info("calculating svd");
    CalculateSingularValueDecomposition();
  }
}

/**
 * @brief load a panel file containing paths of reference coverage files
 * @param pon_fname - path name of reference panel
 * @param on_target - load on-target intervals only if true, else off-targets only
 */
void PanelOfNormals::LoadReferenceCovs(const std::vector<CoverageRecords>& all_covs, bool on_target) {
  std::vector<std::string> all_names;
  std::vector<CoverageRecords> covs_with_target_status;
  size_t pn_intvs = all_covs[0].region.size();
  // filter each normal in the panel by <on_target>-target status
  // also predict the sex of each normal
  for (const auto& covs : all_covs) {
    if (covs.region.size() != pn_intvs) {
      Logging::Error("panel of normals coverage files do not have equal number of targets!");
      throw std::runtime_error("panel of normals coverage files do not have equal number of targets!");
    }
    _xy_ratios.push_back(covs.GetXYRatio().xy_ratio);
    _sex.push_back(covs.PredictSex());
    all_names.push_back(covs.sample_name);
    covs_with_target_status.push_back(covs.GetOnOrOffTargetRegions(on_target));
    pn_intvs = covs.region.size();
  }
  std::vector<std::string> all_intvs = covs_with_target_status[0].region;
  // create matrix of just <on_target>-target regions, rows=intervals, cols=samples
  arma::mat counts(all_intvs.size(), covs_with_target_status.size());
  size_t j = 0;
  counts.each_col([covs_with_target_status, &j](arma::vec& b) {
    b = covs_with_target_status[j].count;
    j++;
  });
  _kept_names = all_names;
  _kept_intvs = all_intvs;
  _counts = counts;
  // create matrix of both on- and off- target> regions, rows=intervals, cols=samples
  // use this matrix to find outliers
  arma::mat full_counts(all_covs[0].region.size(), covs_with_target_status.size());
  size_t i = 0;
  full_counts.each_col([&all_covs, &i](arma::vec& b) {
    b = all_covs[i].count;
    i++;
  });
  arma::uvec sample_filter1 = IdentifySamplesNonExtremeCov(full_counts);
  FilterSamples(sample_filter1);
}

/**
 * FILTERING STRATEGIES
 */

/**
 * @brief Perform successive filtering steps on both the rows (intervals) and
 * columns (samples) of the counts matrix
 */
void PanelOfNormals::FilterCounts() {
  arma::uvec intv_filter1 = IdentifyIntervalsNonLowCov(_counts);
  FilterIntvs(intv_filter1);
  arma::uvec intv_filter2 = IdentifyIntervalsNonZeroCov(_counts);
  FilterIntvs(intv_filter2);
  arma::uvec intv_filter3 = IdentifyIntervalsNonZeroMedian(_counts);
  FilterIntvs(intv_filter3);
}

/**
 * @brief filter intervals given a list of idxs
 * @param idxs - list of idxs
 * @details when filtering rows, also have to filter the list of intervals
 * so that we can keep track of the names of the filtered rows
 */
void PanelOfNormals::FilterIntvs(const arma::uvec& idxs) {
  _counts = _counts.rows(idxs);
  std::vector<std::string> new_intvs;
  for (auto i : idxs) {
    new_intvs.push_back(_kept_intvs[i]);
  }
  _kept_intvs = new_intvs;
}

/**
 * @brief filter samples given a list of idxs
 * @param idxs - list of idxs
 * @details when filtering columns, also have to filter the list of input * filenames * so that we can keep track of the
 * names of the filtered columns. Also have to * filter predicted sex
 */
void PanelOfNormals::FilterSamples(const arma::uvec& idxs) {
  _counts = _counts.cols(idxs);
  std::vector<Sex> new_sexes;
  std::vector<std::string> new_fnames;
  for (auto i : idxs) {
    new_fnames.push_back(_kept_names[i]);
    new_sexes.push_back(_sex[i]);
  }
  _kept_names = new_fnames;
  _sex = new_sexes;
}

/**
 * @brief identify intervals that are non low coverage
 * @details identifies any intervals where the median fractional coverage exceeds the minimum threshold. This threshold
 * is calculated based on the median of the interval median fractional coverages multiplied by kIntvOutlierMin.
 * @return a std::vector of indices indicating intervals with non low coverage
 */
arma::uvec PanelOfNormals::IdentifyIntervalsNonLowCov(const arma::mat& counts) {
  // Calculate fractional coverage. This normalizes for sequencing depth differences between samples allowing us to
  // compare coverage across samples.
  arma::mat frac_counts = counts;
  for (arma::uword col = 0; col < frac_counts.n_cols; col++) {
    frac_counts.col(col) /= arma::sum(frac_counts.col(col));
  }

  // Per interval median fractional coverage
  arma::vec intv_medians = arma::median(frac_counts, 1);

  // Calculate the minimum coverage threshold for an interval to be considered not low coverage
  f64 min_coverage = arma::median(intv_medians) * kIntvOutlierMin;
  arma::uvec not_low_median_cov = intv_medians > min_coverage;
  return arma::find(not_low_median_cov);
}

/**
 * @brief identify samples with non-extreme coverage
 * @return a std::vector of indices indicating samples with non-extreme coverage
 */
arma::uvec PanelOfNormals::IdentifySamplesNonExtremeCov(const arma::mat& counts) {
  arma::vec sample_means = arma::mean(counts, 0).t();
  f64 median_cov = arma::median(sample_means);
  f64 max_cov = kSampleOutlierMaxPct * median_cov;
  f64 min_cov = kSampleOutlierMinPct * median_cov;
  arma::uvec non_outlier_covs = sample_means < max_cov && sample_means > min_cov;
  return arma::find(non_outlier_covs);
}

/**
 * @brief identify intervals where a sufficient proportion of samples have non-zero coverage
 * @return a std::vector of indices indicating intervals where a sufficient proportion of samples have non-zero coverage
 */
arma::uvec PanelOfNormals::IdentifyIntervalsNonZeroCov(const arma::mat& counts) {
  arma::mat filtered_counts = counts;
  arma::vec intv_medians = arma::median(filtered_counts, 1);
  filtered_counts.for_each([](arma::mat::elem_type& v) { v = v > 0 ? 1.0 : 0.0; });
  arma::vec intvs_perc_samples_w_cov = arma::sum(filtered_counts, 1) / static_cast<f64>(filtered_counts.n_cols);
  arma::uvec above_missing_cov = intvs_perc_samples_w_cov > (1 - kIntvOutlierMinZeroPct);
  return arma::find(above_missing_cov);
}

/**
 * @brief identify intervals where median coverage is non-zero
 * @return a std::vector of indices indicating intervals where median coverage is non-zero
 */
arma::uvec PanelOfNormals::IdentifyIntervalsNonZeroMedian(const arma::mat& counts) {
  arma::vec intv_medians = arma::median(counts, 1);
  return arma::find(intv_medians != 0);
}

/**
 * @brief calculate and set the medians for each interval
 * @details after accounting for sex, find the sex-specific median of each interval
 */
void PanelOfNormals::SetIntervalMedians(const arma::mat& counts,
                                        arma::vec& all_medians,
                                        arma::vec& male_medians,
                                        arma::vec& female_medians) {
  // extract autosome regions from counts, calculate interval medians
  all_medians = arma::median(counts, 1);
  // find all allosome intervals
  std::vector<arma::uword> allosome_idxs;
  for (size_t i = 0; i < _kept_intvs.size(); ++i) {
    if (IsInAllosome(_kept_intvs[i])) {
      allosome_idxs.push_back(i);
    }
  }
  /* gather who is male and female */
  std::vector<arma::uword> male_idxs;
  std::vector<arma::uword> female_idxs;
  for (size_t i = 0; i < _sex.size(); ++i) {
    if (_sex[i] == Sex::kMale) {
      male_idxs.push_back(i);
    } else if (_sex[i] == Sex::kFemale) {
      female_idxs.push_back(i);
    }
  }
  // modify allosome interval medians to be sex-specific
  arma::mat male_counts = counts.cols(arma::uvec(male_idxs));
  arma::mat female_counts = counts.cols(arma::uvec(female_idxs));
  if (male_counts.empty()) {
    Logging::Warn("No males found in panel of normals. Male-specific allosome information will be ignored");
  }
  if (female_counts.empty()) {
    Logging::Warn("No females found in panel of normals. Female-specific allosome information will be ignored");
  }
  male_medians = all_medians;
  female_medians = all_medians;
  for (auto i : allosome_idxs) {
    if (!male_counts.empty()) {
      male_medians[i] = arma::median(male_counts.row(i).t());
    }
    if (!female_counts.empty()) {
      female_medians[i] = arma::median(female_counts.row(i).t());
    }
  }
}

/**
 * @brief standardize count matrix to prepare for SVD
 * @details normalize each count by total sample coverage. Find sex-aware
 * interval medians. Impute any 0-cov elements w/ median cov. Truncate extreme
 * values. log2-normalize and median-center
 */
void PanelOfNormals::StandardizeReferenceCovs() {
  // normalize counts per sample (column)
  _original_counts = _counts;
  arma::mat normed_counts = _counts;
  // fractionalize counts by sample
  normed_counts.each_col([](arma::vec& b) { b = b / arma::sum(b); });
  // get the interval medians
  SetIntervalMedians(
      _original_counts, _original_interval_medians, _original_male_interval_medians, _original_female_interval_medians);
  SetIntervalMedians(
      normed_counts, _normed_interval_medians, _normed_male_interval_medians, _normed_female_interval_medians);
  size_t i = 0;
  normed_counts.each_col([this, &i](arma::vec& b) {
    if (_sex[i] == Sex::kMale) {
      b = b / _normed_male_interval_medians;
    } else if (_sex[i] == Sex::kFemale) {
      b = b / _normed_female_interval_medians;
    } else {
      b = b / _normed_interval_medians;
    }
    i += 1;
  });
  // impute 0-coverages with median of interval
  ImputeCounts(normed_counts);
  // replace any element below the 0.001 or above the 0.999 quantiles with the 0.001 or 0.999 quantile, resp
  TruncateOutlierCounts(normed_counts, 0.001, 0.999);
  _counts = normed_counts;
}

/**
 * @brief impute 0-counts by replacing with interval median of non-zero elements
 */
void PanelOfNormals::ImputeCounts(arma::mat& normed_counts) {
  arma::vec non_zero_interval_medians = arma::vec(normed_counts.n_rows, arma::fill::zeros);
  for (size_t i = 0; i < normed_counts.n_rows; ++i) {
    non_zero_interval_medians[i] = arma::median(arma::nonzeros(normed_counts.row(i).t()));
  }
  arma::uword row_idx = 0;
  normed_counts.each_row([&non_zero_interval_medians, &row_idx](arma::rowvec& b) {
    b.for_each([&non_zero_interval_medians, &row_idx](arma::rowvec::elem_type& x) {
      x = math::IsCloseToZero(x) ? non_zero_interval_medians[row_idx] : x;
    });
    row_idx++;
  });
}

/**
 * @brief truncate extreme count values by bringing upper extemes down to the `hi` quantile, and lower extremes up to
 * the `lo` quantile.
 */
void PanelOfNormals::TruncateOutlierCounts(arma::mat& normed_counts, f64 lo, f64 hi) {
  arma::vec low_hi_counts = arma::quantile(normed_counts.as_col(), arma::vec({lo, hi}));
  normed_counts.for_each([&low_hi_counts](arma::vec::elem_type& b) {
    if (b > low_hi_counts[1]) {
      b = low_hi_counts[1];
    } else if (b < low_hi_counts[0]) {
      b = low_hi_counts[0];
    }
  });
}

/**
 * @brief print unnormalized counts to file
 */
void PanelOfNormals::WriteOriginalCounts(std::ostream& ofs) const {
  assert(_kept_intvs.size() == _original_counts.n_rows);
  assert(_kept_names.size() == _original_counts.n_cols);
  ofs << "Target";
  for (const std::string& f : _kept_names) {
    ofs << " " << f;
  }
  ofs << std::endl;
  for (size_t i = 0; i < _original_counts.n_rows; ++i) {
    ofs << _kept_intvs[i];
    for (size_t j = 0; j < _original_counts.n_cols; ++j) {
      ofs << " " << _original_counts(i, j);
    }
    ofs << std::endl;
  }
}

/**
 * @brief print counts to file
 */
void PanelOfNormals::WriteStandardizedCounts(std::ostream& ofs) const {
  assert(_kept_intvs.size() == _counts.n_rows);
  assert(_kept_names.size() == _counts.n_cols);
  ofs << "Target";
  for (const std::string& f : _kept_names) {
    ofs << " " << f;
  }
  ofs << std::endl;
  for (size_t i = 0; i < _counts.n_rows; ++i) {
    ofs << _kept_intvs[i];
    for (size_t j = 0; j < _counts.n_cols; ++j) {
      ofs << " " << _counts(i, j);
    }
    ofs << std::endl;
  }
}

/**
 * @brief print normalized interval medians to file
 */
void PanelOfNormals::WriteNormedIntervalMedians(std::ostream& ofs) const {
  assert(_kept_intvs.size() == _counts.n_rows);
  assert(_kept_names.size() == _counts.n_cols);
  ofs << "Target" << " " << "Target Median" << std::endl;
  for (size_t i = 0; i < _normed_interval_medians.n_elem; ++i) {
    ofs << _kept_intvs[i] << " " << _normed_interval_medians[i] << std::endl;
  }
}

/**
 * @brief print original count interval medians to file
 */
void PanelOfNormals::WriteOriginalIntervalMedians(std::ostream& ofs) const {
  assert(_kept_intvs.size() == _counts.n_rows);
  assert(_kept_names.size() == _counts.n_cols);
  ofs << "Target" << " " << "Target Median" << std::endl;
  for (size_t i = 0; i < _original_interval_medians.n_elem; ++i) {
    ofs << _kept_intvs[i] << " " << _original_interval_medians[i] << std::endl;
  }
}

/**
 * @brief print predicted sexes of samples to file
 */
void PanelOfNormals::WriteSex(std::ostream& ofs) const {
  for (size_t i = 0; i < _kept_names.size(); ++i) {
    ofs << _kept_names[i] << "," << kSexToChar.at(_sex[i]) << "," << _xy_ratios[i] << std::endl;
  }
}

std::vector<std::string> PanelOfNormals::GetIntervals() const {
  return _kept_intvs;
}

arma::vec PanelOfNormals::GetNormedMaleIntervalMedians() const {
  return _normed_male_interval_medians;
}

arma::vec PanelOfNormals::GetNormedFemaleIntervalMedians() const {
  return _normed_female_interval_medians;
}

arma::vec PanelOfNormals::GetNormedAllIntervalMedians() const {
  return _normed_interval_medians;
}

arma::vec PanelOfNormals::GetNormedSexIntervalMedians(Sex sex) const {
  switch (sex) {
    case Sex::kFemale:
      return _normed_female_interval_medians;
    case Sex::kMale:
      return _normed_male_interval_medians;
    default:
      return _normed_interval_medians;
  }
}

arma::vec PanelOfNormals::GetOriginalSexIntervalMedians(Sex sex) const {
  switch (sex) {
    case Sex::kFemale:
      return _original_female_interval_medians;
    case Sex::kMale:
      return _original_male_interval_medians;
    default:
      return _original_interval_medians;
  }
}

const SingularValueDecomposition& PanelOfNormals::GetSingularValueDecomposition() const {
  return _svd;
}

void PanelOfNormals::CalculateSingularValueDecomposition() {
  arma::mat log2_median_centered_counts = _counts;
  log2_median_centered_counts.each_col([](arma::vec& b) { b = arma::log2(b / arma::median(b)); });
  _svd = SingularValueDecomposition(log2_median_centered_counts);
  _svd.TruncateU(_kept_names.size());
}

const f64 kDenoiseImputeVal = 1e-9;

Observations PanelOfNormals::LogRatio(const CoverageRecords& sample, Sex sex) const {
  // filter sample by normaldb intervals
  Observations log_ratios;
  log_ratios.regions = GetIntervals();
  log_ratios.contigs.resize(log_ratios.regions.size());
  log_ratios.starts.resize(log_ratios.regions.size());
  log_ratios.ends.resize(log_ratios.regions.size());
  size_t i = 0;
  for (const auto& region : GetIntervals()) {
    const auto& [contig, start, end] = ParseRegionString(region);
    log_ratios.contigs[i] = contig;
    log_ratios.starts[i] = start;
    log_ratios.ends[i] = end;
    i++;
  }
  log_ratios.obvs = sample.FilterRegion(log_ratios.regions).count;
  // normalize sample counts
  log_ratios.obvs = log_ratios.obvs / arma::sum(log_ratios.obvs);
  // divide by PoN interval medians (this is the ratio)
  arma::vec medians = GetNormedSexIntervalMedians(sex);
  log_ratios.obvs = log_ratios.obvs / medians;
  // impute missing values
  log_ratios.obvs.for_each([](arma::vec::elem_type& x) { x = x > 0 ? x : kDenoiseImputeVal; });
  // log2 normalize and median-center
  log_ratios.obvs = arma::log2(log_ratios.obvs / arma::median(log_ratios.obvs));
  return log_ratios;
}

// TODO: do we have to use a difference svd for each sex?
Observations PanelOfNormals::DenoiseLogR(const Observations& log_ratios, arma::uword num_eigen) const {
  if (num_eigen > _svd.u.n_cols) {
    num_eigen = _svd.u.n_cols;
  }
  arma::mat projection = _svd.u.head_cols(num_eigen);
  // to denoise, subtract projection from ratios
  Observations ret = log_ratios;
  ret.obvs = log_ratios.obvs - (log_ratios.obvs.t() * projection * projection.t()).t();
  return ret;
}

/**
 * @brief use a sample to provide a filter on what intervals in the reference panel should actually be considered
 * @param sample CoverageRecords of sample
 * @param sex  sex of sample
 * @return list of intervals/regions to keep
 */
std::vector<std::string> PanelOfNormals::FilterIntervalsFromSample(const CoverageRecords& sample,
                                                                   Sex sex,
                                                                   size_t min_target_length,
                                                                   f64 min_panel_median_cov,
                                                                   f64 min_panel_median_and_tumor_cov) const {
  std::vector<std::string> intvs = GetIntervals();
  // make sure we only take the intersection of tumor intervals and ref panel intervals
  CoverageRecords sub_cov = sample.FilterRegion(GetIntervals());
  auto medians = GetOriginalSexIntervalMedians(sex);
  std::vector<std::string> ret;
  for (size_t i = 0; i < intvs.size(); ++i) {
    if (intvs[i] != sub_cov.region[i]) {
      throw std::runtime_error("interval mismatch between Reference Panel and sample");
    }
    const auto& [chrom, start, end] = ParseRegionString(intvs[i]);
    // keep interval if target exceeds some minimun length
    bool interval_large_enough = end - start >= min_target_length;
    auto sample_count = sub_cov.count[i];
    auto panel_count = medians[i];
    // keep interval if the reference panel has enough coverage
    bool sufficient_ref_panel_cov = panel_count >= static_cast<f64>(min_panel_median_cov);
    // keep interval if the combined total of the reference panel and the tumor panel have enough coverage
    bool sufficient_total_counts = panel_count + sample_count >= static_cast<f64>(min_panel_median_and_tumor_cov);
    if (interval_large_enough && sufficient_ref_panel_cov && sufficient_total_counts) {
      ret.push_back(intvs[i]);
    }
  }
  Logging::Info("Filtering removed {} intervals", intvs.size() - ret.size());
  if (ret.empty()) {
    throw std::runtime_error("all intervals were filtered out!");
  }
  return ret;
}

}  // namespace xoos::cnc
