#include "denoise/denoise.h"

#include <armadillo>

#include <xoos/log/logging.h>
#include <xoos/util/file-functions.h>

#include "copy-number-caller/copy-number-caller-options.h"
#include "io/column-names.h"
#include "io/copy-number-caller-default-filenames.h"
#include "panel-of-normals.h"
#include "sex.h"
#include "two-sample-logr/calculate-logrs.h"

namespace xoos::cnc {

/**
 * @brief combines on and off target logrs into a single sorted list
 * @param original_covs
 * @param on_target_logrs
 * @param off_target_logrs
 * @return
 */
ObservationsWOnTarget CombineOnOffTargetLogRs(const std::vector<std::string>& regions,
                                              const Observations& on_target_logrs,
                                              const Observations& off_target_logrs) {
  // generate list of all on and off target intervals
  std::vector<std::string> all_intervals;
  all_intervals.reserve(on_target_logrs.regions.size() + off_target_logrs.regions.size());
  all_intervals.insert(all_intervals.end(), on_target_logrs.regions.begin(), on_target_logrs.regions.end());
  all_intervals.insert(all_intervals.end(), off_target_logrs.regions.begin(), off_target_logrs.regions.end());
  arma::vec all_logrs = arma::join_cols(on_target_logrs.obvs, off_target_logrs.obvs);
  // Create std::vector to store the on- and off-target status of each interval
  std::vector<bool> is_on_target(all_intervals.size());
  std::fill(is_on_target.begin(), is_on_target.begin() + static_cast<s32>(on_target_logrs.regions.size()), true);
  std::fill(is_on_target.begin() + static_cast<s32>(on_target_logrs.regions.size()), is_on_target.end(), false);
  // sort all_intervals by comparing against `regions`
  std::unordered_map<std::string, size_t> region_to_idx;
  for (size_t i = 0; i < regions.size(); ++i) {
    region_to_idx[regions[i]] = i;
  }
  Logging::Info("sorting on- and off-target regions");
  // initialize an array that is just 0->#regions-1
  std::vector<arma::uword> sorted_indexes(all_intervals.size());
  std::iota(sorted_indexes.begin(), sorted_indexes.end(), 0);
  // sort sorted_indexes according to all_intervals' sorting order
  std::stable_sort(
      sorted_indexes.begin(), sorted_indexes.end(), [&region_to_idx, &all_intervals](arma::uword i, arma::uword j) {
        return region_to_idx[all_intervals[i]] < region_to_idx[all_intervals[j]];
      });
  // re-sort the logrs, intervals and on-target statuses
  Observations ret;
  ret.obvs.resize(all_logrs.size());
  ret.regions.resize(all_intervals.size());
  std::vector<bool> ret_is_on_target(is_on_target.size());
  for (size_t i = 0; i < sorted_indexes.size(); ++i) {
    size_t j = sorted_indexes[i];
    ret.obvs[i] = all_logrs[j];
    ret.regions[i] = all_intervals[j];
    ret_is_on_target[i] = is_on_target[j];
  }
  return {ret, ret_is_on_target};
}

void WriteSex(Sex sex, std::ostream& os) {
  os << kSexToChar.at(sex) << std::endl;
}

Observations FilterLogrs(const std::vector<std::string>& intervals_to_keep, Sex sex, const Observations& logrs) {
  Observations ret;
  ret.regions.resize(intervals_to_keep.size());
  ret.contigs.resize(intervals_to_keep.size());
  ret.obvs.resize(intervals_to_keep.size());
  ret.starts.resize(intervals_to_keep.size());
  ret.ends.resize(intervals_to_keep.size());
  size_t i = 0;
  for (const auto& intv : intervals_to_keep) {
    // find the index within logrs containing this kept-interval
    const auto& it = std::find(logrs.regions.begin(), logrs.regions.end(), intv);
    if (it == logrs.regions.end()) {
      throw std::runtime_error("region mismatch between logrs and tumor sample");
    }
    size_t idx = std::distance(logrs.regions.begin(), it);
    // assign appropriate interval and logr to the return object
    ret.regions[i] = intv;
    ret.contigs[i] = logrs.contigs[idx];
    ret.obvs[i] = logrs.obvs[idx];
    ret.starts[i] = logrs.starts[idx];
    ret.ends[i] = logrs.ends[idx];
    i++;
  }
  return ret;
}

DenoiseOut Denoise(const CoverageRecords& tumor_cov,
                   const PanelOfNormals& on_target_ref_pool,
                   const PanelOfNormals& off_target_ref_pool,
                   bool no_filter,
                   size_t min_target_length,
                   f64 min_panel_median_cov,
                   f64 min_panel_median_and_tumor_cov,
                   f64 min_off_target_filter_frac) {
  Sex sex = tumor_cov.PredictSex();
  // on-target processing
  Logging::Info("denoising reference pool on-target intervals");
  Observations on_target_logrs = on_target_ref_pool.LogRatio(tumor_cov, sex);
  Observations on_target_denoised_logrs = on_target_ref_pool.DenoiseLogR(on_target_logrs, 20);
  std::vector<std::string> on_target_intervals_to_keep = on_target_ref_pool.FilterIntervalsFromSample(
      tumor_cov, sex, min_target_length, min_panel_median_cov, min_panel_median_and_tumor_cov);
  if (!no_filter) {
    Logging::Info("filtering on-target logr intervals based on tumor sample counts");
    on_target_logrs = FilterLogrs(on_target_intervals_to_keep, sex, on_target_logrs);
    on_target_denoised_logrs = FilterLogrs(on_target_intervals_to_keep, sex, on_target_denoised_logrs);
  }
  // off-target processing
  Logging::Info("denoising reference pool off-target intervals");
  if (!off_target_ref_pool.empty()) {
    Observations off_target_logrs = off_target_ref_pool.LogRatio(tumor_cov, sex);
    Observations off_target_denoised_logrs = off_target_ref_pool.DenoiseLogR(off_target_logrs, 20);
    bool merge_off_and_on_targets = true;
    // filtering
    if (!no_filter) {
      Logging::Info("filtering off-target logr intervals based on tumor sample counts");
      std::vector<std::string> off_target_intervals_to_keep = off_target_ref_pool.FilterIntervalsFromSample(
          tumor_cov, sex, min_target_length, min_panel_median_cov, min_panel_median_and_tumor_cov);
      // we don't need to include off targets if there are very few of them after filtering
      f64 off_target_filter_frac =
          static_cast<f64>(off_target_intervals_to_keep.size()) /
          static_cast<f64>(on_target_intervals_to_keep.size() + off_target_intervals_to_keep.size());
      if (off_target_filter_frac >= min_off_target_filter_frac) {
        off_target_logrs = FilterLogrs(off_target_intervals_to_keep, sex, off_target_logrs);
        off_target_denoised_logrs = FilterLogrs(off_target_intervals_to_keep, sex, off_target_denoised_logrs);
      } else {
        Logging::Info("not enough valid off-target intervals, ignoring off-targets entirely");
        merge_off_and_on_targets = false;
      }
    }
    if (merge_off_and_on_targets) {
      // merge back on- and off- target logrs
      auto [all_logrs, on_target_status] = CombineOnOffTargetLogRs(tumor_cov.region, on_target_logrs, off_target_logrs);
      auto [all_denoised_logrs, on_target_status2] =
          CombineOnOffTargetLogRs(tumor_cov.region, on_target_denoised_logrs, off_target_denoised_logrs);
      if (on_target_status != on_target_status2) {
        throw std::runtime_error(
            "on-target status mismatch between denoised and original logr observations (something went wrong in "
            "denoise)");
      }
      all_denoised_logrs.PopulateFieldsFromRegions();
      all_logrs.PopulateFieldsFromRegions();
      return DenoiseOut{
          .logrs = all_logrs, .denoised_logrs = all_denoised_logrs, .on_target_status = on_target_status, .sex = sex};
    }
  }
  Logging::Info("not using off-target intervals");
  return DenoiseOut{.logrs = on_target_logrs,
                    .denoised_logrs = on_target_denoised_logrs,
                    .on_target_status = std::vector<bool>(on_target_logrs.regions.size(), true),
                    .sex = sex};
}

void DenoiseMain(const CopyNumberCallerOptions& options) {
  std::ifstream tumor_ifs(options.tumor_coverage_fname.value());
  CoverageRecords tumor_cov(tumor_ifs);
  std::ifstream pon_ifs(options.panel_of_normals_lists);
  std::string line;
  std::vector<CoverageRecords> panel_of_normals;
  while (std::getline(pon_ifs, line)) {
    std::ifstream ifs(line);
    panel_of_normals.emplace_back(ifs);
    panel_of_normals.back().sample_name = line;
  }
  const auto denoised_logrs_out = options.output_dir / kDefaultDenoisedLogRsOutput;
  const auto tumor_sex_out = options.output_dir / kDefaultTumorSexOutput;
  const auto logrs_out = options.output_dir / kDefaultLogRsOutput;
  file::CheckFilePermissionsAndOutputPathExistence({}, {denoised_logrs_out, tumor_sex_out, logrs_out});
  std::ofstream denoised_ofs(denoised_logrs_out);
  if (panel_of_normals.size() == 1) {
    Logging::Info("detected only one normal sample. Cannot denoise");
    FilterAndNormalizeTumorNormal(tumor_cov, panel_of_normals[0]);
    Logging::Info("calculating logR by comparing tumor and single normal sample");
    Observations logrs = CalculateLogrs(tumor_cov, panel_of_normals[0]);
    logrs.Write(denoised_ofs, kColumnLogRatio, false, options.command_line_info);
    std::ofstream sex_out(tumor_sex_out);
    WriteSex(tumor_cov.PredictSex(), sex_out);
  } else {
    PanelOfNormals on_target_ref_pool(panel_of_normals, true);
    PanelOfNormals off_target_ref_pool(panel_of_normals, false);
    DenoiseOut ret = Denoise(tumor_cov,
                             on_target_ref_pool,
                             off_target_ref_pool,
                             options.denoise_options.no_filter,
                             options.denoise_options.min_target_length,
                             options.denoise_options.min_panel_median_cov,
                             options.denoise_options.min_panel_median_and_tumor_cov,
                             options.denoise_options.min_off_target_filter_frac);
    ret.denoised_logrs.Write(denoised_ofs, kColumnLogRatio, false, options.command_line_info);
    if (!logrs_out.empty()) {
      std::ofstream logr_ofs(logrs_out);
      ret.logrs.Write(logr_ofs, kColumnLogRatio, false, options.command_line_info);
    }
    std::ofstream sex_out(tumor_sex_out);
    WriteSex(ret.sex, sex_out);
  }
}
}  // namespace xoos::cnc
