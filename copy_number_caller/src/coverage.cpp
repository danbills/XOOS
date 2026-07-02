#include "coverage.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include <csv.hpp>

#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "io/column-names.h"
#include "utility/utility-functions.h"

namespace xoos::cnc {

/**
 * @brief load a coverage BED file.
 * Expected BED format: Contig\tStart\tEnd\tCounts\tTotalCoverage\tOnTarget (0-based half-open).
 * `#`-prefixed metadata/comments and a `#`-prefixed column header are supported and skipped.
 * @param ifs input stream for coverage BED file
 */
void CoverageRecords::LoadCoverageFile(std::istream& ifs) {
  std::string line;
  std::vector<f64> total_coverages;
  std::vector<f64> counts;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::string contig;
    std::string start_str;
    std::string end_str;
    std::string total_cov_str;
    std::string count_str;
    std::string on_target_str;
    std::istringstream ss(line);
    if (!std::getline(ss, contig, '\t') || !std::getline(ss, start_str, '\t') || !std::getline(ss, end_str, '\t') ||
        !std::getline(ss, count_str, '\t') || !std::getline(ss, total_cov_str, '\t') ||
        !std::getline(ss, on_target_str, '\t')) {
      throw std::runtime_error(
          "Coverage BED file must have 6 columns (Contig, Start, End, Counts, TotalCoverage, OnTarget): \"" + line +
          "\"");
    }
    // reconstruct 1-based region string for internal use
    arma::uword start, end;
    try {
      start = std::stoul(start_str);
      end = std::stoul(end_str);
    } catch (const std::exception& e) {
      throw std::runtime_error("Coverage BED file parse error: " + std::string(e.what()) + " (line: \"" + line + "\")");
    }
    std::string r = contig + ":" + std::to_string(start + 1) + "-" + std::to_string(end);
    f64 t = (total_cov_str == "nan" || total_cov_str == "NA") ? 0 : std::stod(total_cov_str);
    if (t < 0) {
      Logging::Error("Invalid total_coverage value {} for region {}", total_cov_str, r);
      throw std::runtime_error("Invalid total_coverage value in coverage file");
    }
    f64 c = (count_str == "nan" || count_str == "NA") ? 0 : std::stod(count_str);
    if (c < 0) {
      Logging::Error("Invalid counts value {} for region {}", count_str, r);
      throw std::runtime_error("Invalid counts value in coverage file");
    }
    region.push_back(r);
    total_coverages.push_back(t);
    counts.push_back(c);
    on_target.push_back(IsTrueString(on_target_str));
  }
  total_coverage = arma::vec(total_coverages);
  count = arma::vec(counts);
  CalculateAverageCoverages();
  UpdateRegionMap();
}

/**
 * @brief repopulates `average_coverage`
 * @details repopulates the `average_coverage` member std::vector by dividing each
 * member in `total_coverage` by the genomic span derived from the corresponding
 * `region` element
 */
void CoverageRecords::CalculateAverageCoverages() {
  average_coverage = arma::vec(total_coverage.n_elem);
  for (size_t i = 0; i < region.size(); ++i) {
    auto reg = region[i];
    auto [contig, start, end] = ParseRegionString(reg);
    average_coverage[i] = static_cast<f64>(total_coverage[i]) / (static_cast<f64>(end) - static_cast<f64>(start));
  }
}

/**
 * @brief returns new CoverageRecords from a subset of row indexes
 * @param idxs - list of row indexes to subset
 * @return subsetted CoverageRecords
 */
CoverageRecords CoverageRecords::FilterRow(const arma::uvec& idxs) const {
  CoverageRecords new_covs;
  for (auto i : idxs) {
    new_covs.region.push_back(region[i]);
    new_covs.on_target.push_back(on_target[i]);
  }
  new_covs.total_coverage = total_coverage.elem(idxs);
  new_covs.average_coverage = average_coverage.elem(idxs);
  new_covs.count = count.elem(idxs);
  if (!mean_mapping_quality.empty()) {
    new_covs.mean_mapping_quality = mean_mapping_quality.elem(idxs);
  }
  new_covs.UpdateRegionMap();
  return new_covs;
}

/**
 * @brief returns new CoverageRecords from a subset containing given chromosome
 * @param chrom - chromsome to subset on
 * @return subsetted CoverageRecords
 */
CoverageRecords CoverageRecords::FilterChrom(const std::string& chrom) const {
  std::vector<arma::uword> idxs;
  for (size_t i = 0; i < region.size(); i++) {
    if (IsEqualContig(region[i], chrom)) {
      idxs.push_back(i);
    }
  }
  return FilterRow(arma::uvec(idxs));
}

/**
 * @brief returns new CoverageRecords from a subset containing the exact regions specified in regions_to_keep (this
 * function performs an lookup of the region string to retrieve the filtered subset; it does NOT try to look up based on
 * coordinate!
 * @param regions_to_keep - regions to subset on
 * @return subsetted CoverageRecords
 */
CoverageRecords CoverageRecords::FilterRegion(const std::vector<std::string>& regions_to_keep) const {
  std::vector<arma::uword> idxs;
  for (const std::string& region : regions_to_keep) {
    auto it = region_to_row.find(region);
    if (it == region_to_row.end()) {
      throw(std::runtime_error("Could not find " + region + " in CoverageRecords"));
    } else {
      idxs.push_back(it->second);
    }
  }
  return FilterRow(idxs);
}

/**
 * @brief returns new CoverageRecords from on- or off-target regions
 * @param o - if true, then on-targets only. else off-targets.
 * @return subsetted CoverageRecords
 */
CoverageRecords CoverageRecords::GetOnOrOffTargetRegions(bool o) const {
  std::vector<arma::uword> idxs;
  for (size_t i = 0; i < on_target.size(); ++i) {
    if ((o && on_target[i]) || (!o && !on_target[i])) {
      idxs.push_back(i);
    }
  }
  return FilterRow(idxs);
}

/**
 * @brief returns new CoverageRecords from regions containing allosomes
 * @return subsetted CoverageRecords
 */
CoverageRecords CoverageRecords::GetAllosomes() const {
  std::vector<arma::uword> idxs;
  for (size_t i = 0; i < region.size(); ++i) {
    if (IsInAllosome(region[i])) {
      idxs.push_back(i);
    }
  }
  return FilterRow(idxs);
}

/**
 * @brief returns new CoverageRecords from regions containing autosomes
 * @return subsetted CoverageRecords
 */
CoverageRecords CoverageRecords::GetAutosomes() const {
  std::vector<arma::uword> idxs;
  for (size_t i = 0; i < region.size(); ++i) {
    if (!IsInAllosome(region[i])) {
      idxs.push_back(i);
    }
  }
  return FilterRow(idxs);
}

/**
 * @brief print coverarge to file
 */
void CoverageRecords::Write(std::ofstream& ofs, const io::CommandLineInfo& command_line_info) const {
  auto writer = csv::make_tsv_writer(ofs);
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info);
  }
  vec<std::string> column_names{
      kColumnContigAsComment, kColumnStart, kColumnEnd, kColumnCounts, kColumnTotalCoverage, kColumnOnTarget};
  writer << column_names;
  for (size_t i = 0; i < region.size(); ++i) {
    const auto& [contig, start, end] = ParseRegionString(region[i]);
    std::string total_cov_str = (!total_coverage.empty()) ? std::to_string(total_coverage[i]) : "NA";
    std::string count_str = (!count.empty()) ? std::to_string(count[i]) : "NA";
    std::string on_target_str = (!on_target.empty()) ? std::to_string(on_target[i]) : "NA";
    writer << std::make_tuple(contig, start, end, count_str, total_cov_str, on_target_str);
  }
}

// mininum number of elements needed in std::vector to remove outliers
// these are implementation details which do not need to be exposed in the header file
const s32 kMinOutlierNelem = 5;
const f64 kIqrMinQuantile = 0.25;
const f64 kIqrMaxQuantile = 0.75;
const f64 kIqrCoeff = 1.5;

/**
 * @brief uses the inter-quantile range (IQR) method to remove outliers. anything outside the IQR (i.e. below the 25%
 * percentile or above the 75% percentile) will be deemed an outlier. If the span of the IQR times IQR_COEFF is less
 * than MIN_H, then the IQR range is too small and outliers will not be removed.
 * @param v std::vector of doubles from which to remove outliers
 * @return arma::vec v with outliers removed
 */
static arma::vec RemoveOutliers(const arma::vec& v) {
  if (v.n_elem < kMinOutlierNelem) {
    return v;
  }
  arma::vec quants = arma::quantile(v, arma::vec({kIqrMinQuantile, kIqrMaxQuantile}));
  f64 iqr = quants[1] - quants[0];
  f64 h = kIqrCoeff * iqr;
  arma::uvec non_outliers = arma::find(v >= quants[0] - h && v <= quants[1] + h);
  return v.elem(non_outliers);
}

/**
 * @brief get chrX / chrY ratio. helpful for sex determination.
 * Supports both hg38 (chrX/chrY) and hg19 (X/Y) contig names.
 */
XYRatioResult CoverageRecords::GetXYRatio() const {
  CoverageRecords on_target_covs = GetOnOrOffTargetRegions(true);
  arma::vec full_chr_x_cov = on_target_covs.FilterChrom("chrX").average_coverage;
  arma::vec full_chr_y_cov = on_target_covs.FilterChrom("chrY").average_coverage;
  // Fall back to unprefixed contig names (hg19/GRCh37)
  if (full_chr_x_cov.empty()) {
    full_chr_x_cov = on_target_covs.FilterChrom("X").average_coverage;
  }
  if (full_chr_y_cov.empty()) {
    full_chr_y_cov = on_target_covs.FilterChrom("Y").average_coverage;
  }
  arma::vec chr_x_cov = RemoveOutliers(full_chr_x_cov);
  arma::vec chr_y_cov = RemoveOutliers(full_chr_y_cov);
  f64 chr_x_avg_cov;
  f64 chr_y_avg_cov;
  if (chr_x_cov.empty() && chr_y_cov.empty()) {
    return {NAN, NAN, NAN};
  }
  if (chr_y_cov.n_elem < 10) {
    chr_x_avg_cov = arma::median(chr_x_cov);
    chr_y_avg_cov = chr_y_cov.n_elem ? arma::median(chr_y_cov) : 0;
  } else {
    chr_x_avg_cov = arma::mean(chr_x_cov);
    chr_y_avg_cov = chr_y_cov.n_elem ? arma::mean(chr_y_cov) : 0;
  }
  f64 xy_ratio = chr_x_avg_cov / (chr_y_avg_cov + 0.0001);
  return {xy_ratio, chr_x_avg_cov, chr_y_avg_cov};
}

/**
 * @brief Predict sex of a sample through ratio of chrX coverage over chrY coverage.
 * If the ratio > MIN_RATIO_FOR_SEX_DET, then return 'F'(Female)
 * if the ratio is between MIN_RATIO_FOR_SEX_DETR and MIN_RATIO_FOR_SEC_NA, then return 'N' (NA)
 * else return 'M' (Male)
 */
Sex CoverageRecords::PredictSex() const {
  XYRatioResult xy_ratio_result = GetXYRatio();
  f64 xy_ratio = xy_ratio_result.xy_ratio;
  if (std::isnan(xy_ratio)) {
    return Sex::kUnknown;
  }
  if (xy_ratio > kMinRatioForSexDet) {
    return Sex::kFemale;
  } else if (xy_ratio > kMinRatioForSexNA) {
    return Sex::kUnknown;
  } else {
    return Sex::kMale;
  }
}

/**
 * @brief update region-to-row map
 */
void CoverageRecords::UpdateRegionMap() {
  region_to_row.clear();
  for (arma::uword i = 0; i < region.size(); ++i) {
    region_to_row[region[i]] = i;
  }
}

/**
 * @brief sort records by given indices
 */
void CoverageRecords::Sort(const std::vector<arma::uword>& argsort) {
  std::vector<std::string> region_resort(region.size());
  std::vector<bool> on_target_resort(on_target.size());
  for (size_t i = 0; i < region.size(); ++i) {
    region_resort[i] = region[argsort[i]];
    on_target_resort[i] = on_target[argsort[i]];
  }
  region = region_resort;
  on_target = on_target_resort;
  arma::uvec argsort_arma(argsort);
  total_coverage = total_coverage.elem(argsort_arma);
  if (!average_coverage.empty()) {
    average_coverage = average_coverage.elem(argsort_arma);
  }
  count = count.elem(argsort_arma);
  UpdateRegionMap();
}

/**
 * @brief sort records by region
 */
void CoverageRecords::Sort() {
  std::vector<arma::uword> argsort(region.size());
  std::iota(argsort.begin(), argsort.end(), 0);
  std::stable_sort(argsort.begin(), argsort.end(), [this](size_t x, size_t y) { return region[x] < region[y]; });
  Sort(argsort);
}

/**
 * @brief sort records by given regions
 */
void CoverageRecords::Sort(const std::vector<std::string>& ordered_regions) {
  std::vector<arma::uword> argsort(region.size());
  size_t j = 0;
  for (const auto& r : ordered_regions) {
    auto it = region_to_row.find(r);
    if (it == region_to_row.end()) {
      Logging::Debug("{} does not appear in coverage object. skipping...", r);
      continue;
    }
    argsort[j++] = it->second;
  }
  if (j != region.size()) {
    throw std::runtime_error("Baits file is missing a region that was encountered in coverage file");
  }
  Sort(argsort);
}

bool VerifyRegions(const std::vector<std::string>& lhs, const std::vector<std::string>& rhs) {
  // we assume that the baits and the coverage records have completely the same regions (no subsets)
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (size_t i = 0; i < lhs.size(); ++i) {
    if (lhs[i] != rhs[i]) {
      Logging::Error("region mismatch: bait {} vs coverage {} records", lhs[i], rhs[i]);
      return false;
    }
  }
  return true;
}

}  // namespace xoos::cnc
