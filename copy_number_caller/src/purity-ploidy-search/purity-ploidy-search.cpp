#include "purity-ploidy-search/purity-ploidy-search.h"

#include <cmath>
#include <stdexcept>

#include <csv.hpp>

#include <xoos/enum/enum-util.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>
#include <xoos/stats/copy-number-stats.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>

#include "copy-number-caller/copy-number-caller-options.h"
#include "io/copy-number-caller-default-filenames.h"
#include "likelihood/baf-tools.h"
#include "likelihood/likelihood.h"
#include "segmentation/interval-trees.h"
#include "segmentation/read-segments.h"

namespace xoos::cnc {
const size_t kLikelihoodMaxCN1 = 7;  // for CN likelihood estimation in purity-ploidy search
const f64 kLikelihoodRoStep = 0.01;  // Step size for purity
const f64 kLikelihoodRoMin = 0.15;   // Minimum purity to consider
const f64 kLikelihoodRoMax = 0.95;   // Maximum purity to consider
const size_t kLikelihoodRoNumSteps = static_cast<size_t>((kLikelihoodRoMax - kLikelihoodRoMin) / kLikelihoodRoStep) + 1;
const f64 kLikelihoodPhiStep = 0.1;  // Step size for ploidy
const f64 kLikelihoodPhiMin = 1.4;   // Minimum ploidy to consider
const f64 kLikelihoodPhiMax = 6;     // Maximum ploidy to consider
const size_t kLikelihoodPhiNumSteps =
    static_cast<size_t>((kLikelihoodPhiMax - kLikelihoodPhiMin) / kLikelihoodPhiStep) + 1;

const f64 kMinIntervalDensity = 0.5;
constexpr s32 kPercentMultiplier = 100;
const size_t kMinCandidateHetDelSegments = 5;
const f64 kHetDelExpectedTcnMin = 0.9;
const f64 kHetDelExpectedTcnMax = 1.1;
const f64 kHetDelMaxMBAF = 0.45;
const f64 kWGDMaxPeakDifference = 0.1;
const f64 kMinPurityToEvaluate = 0.3;
const size_t kNumTopPeaks = 2;

/**
 * @brief Use global tumor purity/ploidy estimates (user-provided or optionally estimate those here) to calculate the
 * joint logr-BAF likelihood of each segment. In calculating this joint likelihood we also return the most likely total
 * copy number of the segment as well as the the most likely minor copy number. See docs/likelihood.md for more
 * information
 * @param params
 */
void PurityPloidySearchMain(const CopyNumberCallerOptions& options) {
  std::ifstream logr_ifstream(options.logrs_fname.value());
  Observations logrs = ReadObservations(logr_ifstream);

  // create output files
  const auto purity_ploidy_out_fname = options.output_dir / kDefaultPurityPloidyOutput;
  const auto purity_ploidy_grid_out_fname = options.output_dir / kDefaultPurityPloidyGridOutput;
  file::CheckFilePermissionsAndOutputPathExistence({}, {purity_ploidy_out_fname, purity_ploidy_grid_out_fname});

  std::optional<Observations> ref_ads, alt_ads;
  segmentation::SegmentType segment_type;
  if (!options.bafs_fname.has_value()) {
    segment_type = SegmentType::kLogROnly;
  } else {
    auto ref_alts = GetDepthsFromFile(options.bafs_fname.value());
    ref_ads = std::move(ref_alts.ref_obvs);
    alt_ads = std::move(ref_alts.alt_obvs);
    segment_type = SegmentType::kBaf;
  }
  auto segments = segmentation::ReadSegments(options.segments_fname.value(), segment_type);
  // Recompute mean_logr and num_obs from raw per-locus observations so that all downstream
  // steps (ShouldSkipSegment, EvaluatePurityPloidySolution, grid search) use consistent
  // full-precision values instead of the rounded MeanLogRatio from the SEG file.
  segmentation::PopulateGenomicSegmentOptionalFields(segments, logrs);
  const PurityPloidySearchOut purity_ploidy_out =
      PurityPloidySearchWithWGDQC(segments,
                                  logrs,
                                  ref_ads,
                                  alt_ads,
                                  options.purity_ploidy_search_options,
                                  purity_ploidy_grid_out_fname,
                                  options.sample_metadata_options.sex.value(),
                                  options.command_line_info);
  WritePurityPloidy(purity_ploidy_out.purity,
                    purity_ploidy_out.ploidy,
                    purity_ploidy_out_fname,
                    options.command_line_info,
                    purity_ploidy_out.total_copy_number_prior);
  if (ref_ads.has_value() && alt_ads.has_value()) {
    CheckExtremeBafProportionAndWarn(ref_ads.value(), alt_ads.value(), purity_ploidy_out.purity);
  }
}

PurityPloidySearchOut PurityPloidySearch(const std::vector<segmentation::GenomicSegment>& segments,
                                         const Observations& logrs,
                                         const std::optional<Observations>& ref_obvs,
                                         const std::optional<Observations>& alt_obvs,
                                         const f64 seg_max_logr,
                                         const size_t seg_min_num_logrs,
                                         const size_t seg_min_num_snps,
                                         const fs::path& grid_out_fname,
                                         const Sex sex,
                                         const io::CommandLineInfo& command_line_info,
                                         const bool use_diploid_bias_prior) {
  f64 purity, ploidy;
  arma::mat ppgrid_joint, ppgrid_logr, ppgrid_baf;
  const TotalCopyNumberPrior prior =
      use_diploid_bias_prior ? TotalCopyNumberPrior::kDiploidBias : TotalCopyNumberPrior::kUniform;
  if (ref_obvs.has_value() && alt_obvs.has_value()) {
    std::tie(ppgrid_joint, ppgrid_logr, ppgrid_baf) = GetPurityPloidyGridSerialLogRBAF(segments,
                                                                                       logrs,
                                                                                       ref_obvs.value(),
                                                                                       alt_obvs.value(),
                                                                                       seg_max_logr,
                                                                                       seg_min_num_logrs,
                                                                                       seg_min_num_snps,
                                                                                       sex,
                                                                                       use_diploid_bias_prior);
    std::tie(purity, ploidy) = EstimatePurityPloidyByTotalMaximum(ppgrid_joint);
    WritePurityPloidyGrid(ppgrid_joint, ppgrid_logr, ppgrid_baf, grid_out_fname, command_line_info, prior);
  } else {
    ppgrid_joint =
        GetPurityPloidyGridLogROnly(segments, logrs, seg_max_logr, seg_min_num_logrs, sex, use_diploid_bias_prior);
    std::tie(purity, ploidy) = EstimatePurityPloidyByTotalMaximum(ppgrid_joint);
    WritePurityPloidyGrid(ppgrid_joint, {}, {}, grid_out_fname, command_line_info, prior);
  }
  return {.purity = purity,
          .ploidy = ploidy,
          .joint_ll_grid = ppgrid_joint,
          .logr_ll_grid = ppgrid_logr,
          .baf_ll_grid = ppgrid_baf,
          .total_copy_number_prior = prior};
}

/**
 * @brief Checks if a segment should be skipped based on mean_logr, num_obs, and optionally num_snps.
 * @param seg The segment to check.
 * @param seg_max_logr Maximum allowed mean_logr.
 * @param seg_min_num_logrs Minimum required number of logRs.
 * @param seg_min_num_snps Minimum required number of SNPs.
 * @param check_seg_min_snps Whether to check for minimum required number of SNPs.
 * @return bool True if the segment should be skipped, false otherwise.
 */
bool ShouldSkipSegment(const GenomicSegment& seg,
                       f64 seg_max_logr,
                       size_t seg_min_num_logrs,
                       std::optional<size_t> seg_min_num_snps,
                       bool check_seg_min_snps) {
  if (seg.mean_logr > seg_max_logr) {
    Logging::Info("Skipping segment {}:{}-{} because mean logr value is greater than {}",
                  seg.contig,
                  seg.start,
                  seg.end,
                  seg_max_logr);
    return true;
  }

  if (seg.num_obs.value() < seg_min_num_logrs) {
    Logging::Info(
        "Skipping segment {}:{}-{} because it has < {} logr", seg.contig, seg.start, seg.end, seg_min_num_logrs);
    return true;
  }
  if (check_seg_min_snps) {
    if (seg.num_snps.value() < seg_min_num_snps.value()) {
      Logging::Info("Skipping segment {}:{}-{} because it has < {} SNPs",
                    seg.contig,
                    seg.start,
                    seg.end,
                    seg_min_num_snps.value());
      return true;
    }
  }
  return false;
}

/**
 * @brief Checks whether a segment has low interval density relative to its length.
 * A segment passes (returns false) when the total base coverage of the logR intervals that overlap the segment is at
 * least kMinIntervalDensity * segment_length. Segments that fail this check are large regions with very sparse interval
 * coverage and should be excluded from purity/ploidy estimation.
 * @param logr_idxs Indices into logrs for the intervals that overlap the segment.
 * @param logrs All logR observations (provides genomic positions via starts/ends).
 * @param seg The segment being evaluated.
 * @return true if the interval density is too low (segment should be skipped), false otherwise.
 */
bool HasLowIntervalDensity(const vec<size_t>& logr_idxs, const Observations& logrs, const GenomicSegment& seg) {
  if (seg.end == seg.start) {
    return false;
  }
  size_t total_covered_bases = 0;
  for (const size_t idx : logr_idxs) {
    total_covered_bases += logrs.ends[idx] - logrs.starts[idx];
  }
  const auto segment_length = static_cast<f64>(seg.end - seg.start);
  if (static_cast<f64>(total_covered_bases) < kMinIntervalDensity * segment_length) {
    Logging::Info("Skipping segment {}:{}-{} because interval coverage ({} bp) is below {}% of segment length ({} bp)",
                  seg.contig,
                  seg.start,
                  seg.end,
                  total_covered_bases,
                  static_cast<s32>(kMinIntervalDensity * kPercentMultiplier),
                  static_cast<size_t>(segment_length));
    return true;
  }
  return false;
}

/**
 * @brief Builds a grid of LogR likelihoods over each combination of purity and ploidy, where purity and ploidy steps
 * are pre-defined
 * @param segments
 * @param obvs
 * @return arma::mat - max LogR likelihood over all purity-ploidy steps
 */
arma::mat GetPurityPloidyGridLogROnly(const std::vector<segmentation::GenomicSegment>& segments,
                                      const Observations& obvs,
                                      f64 seg_max_logr,
                                      const size_t seg_min_num_logrs,
                                      const Sex sex,
                                      const bool use_diploid_bias_prior) {
  // iniitialize a 2D grid - purity, ploidy for each segment
  std::vector<arma::mat> grids;
  grids.reserve(segments.size());
  for (size_t i = 0; i < segments.size(); ++i) {
    grids.emplace_back(kLikelihoodRoNumSteps, kLikelihoodPhiNumSteps);
  }
  segmentation::IntervalTrees trees(obvs);
  f64 sd = GetMeanSegmentLogRSD(segments, obvs, trees);
  for (size_t seg_idx = 0; seg_idx < segments.size(); ++seg_idx) {
    // skip sex chromosomes
    const GenomicSegment& seg = segments[seg_idx];
    Logging::Info("building grid for segment #{}", seg_idx);
    if (SkipBasedOnModeAndSex(seg, LikelihoodMode::kGermline, sex)) {
      Logging::Info("skipping segment #{} {}:{}-{} because due to sex of sample and mode",
                    seg_idx,
                    seg.contig,
                    seg.start,
                    seg.end);
      continue;
    }
    bool expect_haploid = ExpectHaploid(seg, LikelihoodMode::kGermline, sex);
    auto seg_obvs_idxs = trees.LookUp(seg.contig, seg.start, seg.end);
    if (seg_obvs_idxs.empty()) {
      Logging::Info("skipping segment #{} {}:{}-{} because it has no observations. Probably a seed segment",
                    seg_idx,
                    seg.contig,
                    seg.start,
                    seg.end);
      continue;
    }
    if (HasLowIntervalDensity(seg_obvs_idxs, obvs, seg)) {
      continue;
    }
    auto skip_segment = ShouldSkipSegment(seg, seg_max_logr, seg_min_num_logrs, std::nullopt, false);
    if (skip_segment) {
      continue;
    }
    auto seg_obvs = obvs.FilterByIdxs(seg_obvs_idxs);
    auto& grid = grids[seg_idx];
    // TODO: if we decide to parallelize, then we have to initialize this per loop
    arma::vec logr_likelihoods(kLikelihoodMaxCN1 + 1);
    logr_likelihoods.fill(kNegativeInf);
    arma::mat::row_col_iterator it = grid.begin_row_col();
    for (; it != grid.end_row_col(); ++it) {
      size_t row_ind = it.row();
      f64 ro = kLikelihoodRoMin + (static_cast<f64>(row_ind) * kLikelihoodRoStep);
      size_t col_ind = it.col();
      f64 phi = kLikelihoodPhiMin + (static_cast<f64>(col_ind) * kLikelihoodPhiStep);
      // calculate the most likely copy number and its likelihood for this cell
      for (size_t i = 0; i < logr_likelihoods.size(); ++i) {
        auto tcn = static_cast<f64>(i);
        auto n_cns = static_cast<f64>(logr_likelihoods.size());
        const f64 weight = (use_diploid_bias_prior && i == 2) ? 2.0 / n_cns : 1.0 / n_cns;
        logr_likelihoods[i] = weight * SegmentSummarizedLogRLikelihood(seg_obvs.obvs, sd, ro, phi, tcn, expect_haploid);
      }
      *it = arma::max(logr_likelihoods);
      logr_likelihoods.fill(0);  // reset the std::vector
    }
  }
  // Initialize a matrix to hold the sum of all grids
  arma::mat sum_of_grids(kLikelihoodRoNumSteps, kLikelihoodPhiNumSteps, arma::fill::zeros);
  // Iterate over the grids std::vector and add each grid to sum_of_grids
  for (const auto& grid : grids) {
    sum_of_grids += grid;
  }
  // sum up all the girds
  return {sum_of_grids};
}

/**
 * @brief Builds a grid of joint LogR-BAF likelihoods over each combination of purity and ploidy, where purity and
 * ploidy steps are pre-defined
 * @param segments
 * @param obvs
 * @return matrices of joint-likelihoods, logrs, bafs for each purity-ploidy combination
 */
std::tuple<arma::mat, arma::mat, arma::mat> GetPurityPloidyGridSerialLogRBAF(
    const std::vector<GenomicSegment>& segments,
    const Observations& logrs,
    const Observations& ref_ads,
    const Observations& alt_ads,
    const f64 seg_max_logr,
    const size_t seg_min_num_logrs,
    const size_t seg_min_num_snps,
    const Sex sex,
    const bool use_diploid_bias_prior) {
  // iniitialize a 2D grid - purity, ploidy for each segment
  std::vector<arma::mat> joint_ll_grids;
  std::vector<arma::mat> logr_ll_grids;
  std::vector<arma::mat> baf_ll_grids;
  joint_ll_grids.reserve(segments.size());
  logr_ll_grids.reserve(segments.size());
  baf_ll_grids.reserve(segments.size());
  for (size_t i = 0; i < segments.size(); ++i) {
    joint_ll_grids.emplace_back(kLikelihoodRoNumSteps, kLikelihoodPhiNumSteps);
    logr_ll_grids.emplace_back(kLikelihoodRoNumSteps, kLikelihoodPhiNumSteps);
    baf_ll_grids.emplace_back(kLikelihoodRoNumSteps, kLikelihoodPhiNumSteps);
  }
  segmentation::IntervalTrees logr_trees(logrs);
  IntervalTrees variant_trees(ref_ads);
  // TODO: if we decide to parallelize, then we have to initialize this per loop
  arma::vec logr_likelihoods(kLikelihoodMaxCN1 + 1);
  arma::vec baf_likelihoods(kLikelihoodMaxCN1 + 1);
  f64 mean_segment_logr_sd = GetMeanSegmentLogRSD(segments, logrs, logr_trees);
  // The mirror process reduces the variances of the BAFs. We f64 the mean segment BAF SD to account for this
  // we also w
  f64 mean_segment_mbaf_sd = 2 * GetMeanSegmentMBAFSD(segments, ref_ads, alt_ads, variant_trees);
  for (size_t seg_idx = 0; seg_idx < segments.size(); ++seg_idx) {
    // skip sex chromosomes
    const GenomicSegment& seg = segments[seg_idx];
    Logging::Info("building grid for segment #{}", seg_idx);
    if (SkipBasedOnModeAndSex(seg, LikelihoodMode::kSomatic, sex)) {
      Logging::Info("skipping segment #{} {}:{}-{} because due to sex of sample and mode",
                    seg_idx,
                    seg.contig,
                    seg.start,
                    seg.end);
      continue;
    }
    bool expect_haploid = ExpectHaploid(seg, LikelihoodMode::kSomatic, sex);
    auto& joint_ll_grid = joint_ll_grids[seg_idx];
    auto& logr_ll_grid = logr_ll_grids[seg_idx];
    auto& baf_ll_grid = baf_ll_grids[seg_idx];
    std::vector<size_t> logr_idxs = logr_trees.LookUp(seg.contig, seg.start, seg.end);
    if (logr_idxs.empty()) {
      Logging::Info(
          "skipping segment #{} {}:{}-{} because it has no logr observations", seg_idx, seg.contig, seg.start, seg.end);
      continue;
    }
    if (HasLowIntervalDensity(logr_idxs, logrs, seg)) {
      continue;
    }
    auto skip_segment = ShouldSkipSegment(seg, seg_max_logr, seg_min_num_logrs, seg_min_num_snps);
    if (skip_segment) {
      continue;
    }

    Observations logr_subset = logrs.FilterByIdxs(logr_idxs);
    std::vector<size_t> variant_idxs = variant_trees.LookUp(seg.contig, seg.start, seg.end);
    Observations ref_ad_subset = ref_ads.FilterByIdxs(variant_idxs);
    Observations alt_ad_subset = alt_ads.FilterByIdxs(variant_idxs);
    arma::mat::const_row_col_iterator it = joint_ll_grid.begin_row_col();
    std::optional<f64> seg_mbaf = std::nullopt;
    if (!ref_ad_subset.obvs.empty()) {
      f64 mbaf = GetPeakOfBAFDistributionThenMirror(ref_ad_subset.obvs, alt_ad_subset.obvs);
      seg_mbaf = std::isnan(mbaf) ? std::nullopt : std::optional<f64>(mbaf);
    }
    for (; it != joint_ll_grid.end_row_col(); ++it) {
      size_t row_ind = it.row();
      f64 purity = kLikelihoodRoMin + (static_cast<f64>(row_ind) * kLikelihoodRoStep);
      size_t col_ind = it.col();
      f64 ploidy = kLikelihoodPhiMin + (static_cast<f64>(col_ind) * kLikelihoodPhiStep);
      GenomicSegment segment;
      segment.mbaf = seg_mbaf;
      segment = SegmentSerialSummarizedLogRBAFLikelihood(segment,
                                                         logr_subset.obvs,
                                                         purity,
                                                         ploidy,
                                                         mean_segment_logr_sd,
                                                         mean_segment_mbaf_sd,
                                                         use_diploid_bias_prior,
                                                         expect_haploid,
                                                         logr_likelihoods,
                                                         baf_likelihoods);
      joint_ll_grid(row_ind, col_ind) = segment.joint_likelihood.value();
      logr_ll_grid(row_ind, col_ind) = segment.logr_likelihood.value();
      baf_ll_grid(row_ind, col_ind) = segment.baf_likelihood.value_or(NAN);
    }
  }
  // Initialize a matrix to hold the sum of all grids
  arma::mat sum_joint_ll_grids(kLikelihoodRoNumSteps, kLikelihoodPhiNumSteps, arma::fill::zeros);
  arma::mat sum_logr_ll_grids(kLikelihoodRoNumSteps, kLikelihoodPhiNumSteps, arma::fill::zeros);
  arma::mat sum_baf_ll_grids(kLikelihoodRoNumSteps, kLikelihoodPhiNumSteps, arma::fill::zeros);
  // Iterate over the grids std::vector and add each grid to sum_of_grids
  for (const auto& joint_ll_grid : joint_ll_grids) {
    sum_joint_ll_grids += joint_ll_grid;
  }
  for (const auto& logr_ll_grid : logr_ll_grids) {
    sum_logr_ll_grids += logr_ll_grid;
  }
  for (auto& baf_ll_grid : baf_ll_grids) {
    baf_ll_grid.replace(arma::datum::nan, 0);
    sum_baf_ll_grids += baf_ll_grid;
  }
  // sum up all the girds
  return {sum_joint_ll_grids, sum_logr_ll_grids, sum_baf_ll_grids};
}

/**
 * @brief given a purity ploidy grid, return the purity ploidy associated with the absolute maximum of the
 * grid
 * @param purity_ploidy_grid
 * @return
 */
tuple<f64, f64> EstimatePurityPloidyByTotalMaximum(const arma::mat& purity_ploidy_grid) {
  size_t idx_max_ll = purity_ploidy_grid.index_max();
  arma::uvec coords_max_ll = arma::ind2sub(arma::size(purity_ploidy_grid), idx_max_ll);
  f64 purity = kLikelihoodRoMin + kLikelihoodRoStep * static_cast<f64>(coords_max_ll(0));
  f64 ploidy = kLikelihoodPhiMin + kLikelihoodPhiStep * static_cast<f64>(coords_max_ll(1));
  return {purity, ploidy};
}

tuple<f64, f64> EstimatePurityPloidyByLocalMaxima(const arma::mat& purity_ploidy_grid,
                                                  const Observations& ref_ads,
                                                  const Observations& alt_ads) {
  // find local maxima with simulated annealing
  throw std::runtime_error("EstimatePurityPloidyByLocalMaxima not yet implemented");
  return {};
}

/**
 * @brief write joint, logr and BAF likelihoods from the purity ploidy grid to disk
 * @param purity_ploidy_grid
 * @param out_fname
 */
void WritePurityPloidyGrid(const arma::mat& purity_ploidy_grid_joint_lls,
                           const std::optional<arma::mat>& purity_ploidy_grid_logr_lls,
                           const std::optional<arma::mat>& purity_ploidy_grid_baf_lls,
                           const fs::path& out_fname,
                           const io::CommandLineInfo& command_line_info,
                           std::optional<TotalCopyNumberPrior> prior) {
  std::ofstream ofs(out_fname);
  if (!ofs.is_open()) {
    Logging::Error("Failed to open purity ploidy grid output file: {}", out_fname.string());
    throw std::runtime_error("Failed to open purity ploidy grid output file: " + out_fname.string());
  }
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info);
  }
  if (prior.has_value()) {
    ofs << "#TotalCopyNumberPrior=" << enum_util::FormatEnumName(prior.value()) << "\n";
  }
  ofs << kPurityPloidyGridOutPurity << "\t" << kPurityPloidyGridOutPloidy << "\t"
      << kPurityPloidyGridOutTotalJointLikelihood << "\t" << kPurityPloidyGridOutTotalLogRatioLikelihood << "\t"
      << kPurityPloidyGridOutTotalBAlleleFrequencyLikelihood << std::endl;
  f64 purity_val;
  f64 ploidy_val;
  for (arma::uword row_idx = 0; row_idx < purity_ploidy_grid_joint_lls.n_rows; ++row_idx) {
    // Determine corresponding purity value of the row index
    purity_val = kLikelihoodRoMin + (static_cast<f64>(row_idx) * kLikelihoodRoStep);
    for (arma::uword col_idx = 0; col_idx < purity_ploidy_grid_joint_lls.n_cols; ++col_idx) {
      // Determine corresponding ploidy value of the column index
      ploidy_val = kLikelihoodPhiMin + static_cast<f64>(col_idx) * kLikelihoodPhiStep;
      ofs << purity_val << "\t" << ploidy_val << "\t" << purity_ploidy_grid_joint_lls(row_idx, col_idx);
      if (purity_ploidy_grid_logr_lls.has_value()) {
        ofs << "\t" << purity_ploidy_grid_logr_lls.value()(row_idx, col_idx);
      } else {
        ofs << "\t" << "nan";
      }
      if (purity_ploidy_grid_baf_lls.has_value()) {
        ofs << "\t" << purity_ploidy_grid_baf_lls.value()(row_idx, col_idx);
      } else {
        ofs << "\t" << "nan";
      }
      ofs << std::endl;
    }
  }
}

/**
 * @brief write purity and ploidy to disk
 * @param purity
 * @param ploidy
 * @param out_fname
 */
void WritePurityPloidy(const f64 purity,
                       const f64 ploidy,
                       const fs::path& out_fname,
                       const io::CommandLineInfo& command_line_info,
                       std::optional<TotalCopyNumberPrior> prior) {
  std::ofstream ofs(out_fname);
  if (!ofs.is_open()) {
    Logging::Error("Failed to open purity ploidy output file: {}", out_fname.string());
    throw std::runtime_error("Failed to open purity ploidy output file: " + out_fname.string());
  }
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info);
  }
  if (prior.has_value()) {
    ofs << "#TotalCopyNumberPrior=" << enum_util::FormatEnumName(prior.value()) << "\n";
  }
  ofs << kPurityPloidyOutPurity << "\t" << kPurityPloidyOutPloidy << std::endl;
  ofs << purity << "\t" << ploidy << "\n";
}

namespace {

/**
 * @brief Identifies candidate clonal heterozygous deletion segments for WGD QC.
 *
 * Filters segments by: autosome membership, minimum logR/SNP observations, interval density,
 * and expected total copy number in the het-del range.
 */
vec<std::pair<const GenomicSegment*, f64>> FindCandidateHetDelSegments(const vec<GenomicSegment>& segments,
                                                                       const Observations& logrs,
                                                                       const RefAltObservations& ref_alt_depths,
                                                                       const f64 purity,
                                                                       const f64 ploidy,
                                                                       const PurityPloidySearchOptions& options) {
  const segmentation::IntervalTrees logr_trees(logrs);
  const segmentation::IntervalTrees variant_trees(ref_alt_depths.ref_obvs);
  vec<std::pair<const GenomicSegment*, f64>> candidates;
  for (const auto& seg : segments) {
    if (IsInAllosome(seg.contig) || !seg.mean_logr.has_value()) {
      continue;
    }
    if (!seg.num_obs.has_value() || seg.num_obs.value() < options.seg_min_num_logrs) {
      continue;
    }
    if (!seg.num_snps.has_value() || seg.num_snps.value() < options.seg_min_num_snps) {
      continue;
    }
    const std::vector<size_t> logr_idxs = logr_trees.LookUp(seg.contig, seg.start, seg.end);
    if (HasLowIntervalDensity(logr_idxs, logrs, seg)) {
      continue;
    }
    const f64 expected_tcn = stats::ExpectedTotalCopyNumber(purity, ploidy, seg.mean_logr.value(), false);
    if (!(expected_tcn >= kHetDelExpectedTcnMin && expected_tcn <= kHetDelExpectedTcnMax)) {
      continue;
    }
    const vec<size_t> variant_idxs = variant_trees.LookUp(seg.contig, seg.start, seg.end);
    if (variant_idxs.empty()) {
      continue;
    }
    const Observations ref_subset = ref_alt_depths.ref_obvs.FilterByIdxs(variant_idxs);
    const Observations alt_subset = ref_alt_depths.alt_obvs.FilterByIdxs(variant_idxs);
    const f64 seg_mbaf = GetPeakOfBAFDistributionThenMirror(ref_subset.obvs, alt_subset.obvs);
    if (std::isnan(seg_mbaf)) {
      continue;
    }
    if (seg_mbaf >= kHetDelMaxMBAF) {
      Logging::Info("WGD QC: excluding candidate segment {}:{}-{} because mBAF {:.4f} > {:.2f}",
                    seg.contig,
                    seg.start,
                    seg.end,
                    seg_mbaf,
                    kHetDelMaxMBAF);
      continue;
    }
    candidates.emplace_back(&seg, seg_mbaf);
  }
  return candidates;
}

}  // namespace

/**
 * @brief Checks whether purity/ploidy solution is valid using a diploid bias prior
 *
 * The algorithm may converge on ploidy~2 even though the real ploidy is ~4 (whole genome doubling, WGD) when using a
 * diploid bias prior for total copy number. This function detects that scenario by inspecting the mirrored BAF
 * distribution of putative clonal heterozygous-deletion segments (expected total copy number ≈ 1).
 *
 * In a true WGD tumour those segments are actually copy-neutral LOH (TCN≈2), so their mBAF values cluster near 0 and
 * near 0.5 (assuming 100% purity) simultaneously, producing two distinct peaks in the KDE. A non-WGD tumour produces a
 * single peak for these deletions. If two peaks separated by more than @p kWGDMaxPeakDifference are detected the
 * function returns TotalCopyNumberPrior::kUniform, signalling that the search should be re-run with a uniform
 * copy-number prior.
 *
 * Early-exit conditions that accept the diploid solution without further testing:
 *   - purity < kMinPurityToEvaluate
 *   - fewer than kMinCandidateHetDelSegments segments pass all filters
 *   - fewer than kMinCandidateHetDelSegments of those segments yield a valid mBAF estimate
 *   - the KDE of candidate-segment mBAFs has fewer than two local maxima
 * @param segments All genomic segments for the sample.
 * @param ref_alt_depths Per-locus reference and alternate allele depths (used to look up BAF values).
 * @param purity Estimated tumor purity from the purity/ploidy search.
 * @param ploidy Estimated tumor ploidy from the purity/ploidy search.
 * @param options Purity/ploidy search options (segment filter thresholds).
 * @return TotalCopyNumberPrior::kDiploidBias if the diploid solution is accepted,
 *         TotalCopyNumberPrior::kUniform if it is rejected (WGD suspected).
 */
TotalCopyNumberPrior EvaluatePurityPloidySolution(
    const std::vector<GenomicSegment>&
        segments,  // NOSONAR(S831): intentionally internal, tested via forward declaration
    const Observations& logrs,
    const RefAltObservations& ref_alt_depths,
    const f64 purity,
    const f64 ploidy,
    const PurityPloidySearchOptions& options) {
  if (purity < kMinPurityToEvaluate) {
    Logging::Info(
        "WGD QC: purity {:.2f} is below threshold {:.2f}, accepting diploid solution", purity, kMinPurityToEvaluate);
    return TotalCopyNumberPrior::kDiploidBias;
  }

  const auto candidate_segments = FindCandidateHetDelSegments(segments, logrs, ref_alt_depths, purity, ploidy, options);

  Logging::Info("WGD QC: {} candidate clonal het-del segments identified (need {})",
                candidate_segments.size(),
                kMinCandidateHetDelSegments);
  for (const auto& [seg, mbaf] : candidate_segments) {
    Logging::Info(
        "WGD QC:   candidate segment {}:{}-{} num_obs={} num_snps={} mean_logr={:.6f} mean_dh={:.4f} mbaf={:.4f}",
        seg->contig,
        seg->start,
        seg->end,
        seg->num_obs.value_or(0),
        seg->num_snps.value_or(0),
        seg->mean_logr.value(),
        seg->mean_dh.value_or(NAN),
        mbaf);
  }

  if (candidate_segments.size() < kMinCandidateHetDelSegments) {
    Logging::Info(
        "WGD QC: only {} candidate clonal het-del segments (need {}), accepting purity/ploidy solution "
        "(purity={:.2f}, ploidy={:.2f})",
        candidate_segments.size(),
        kMinCandidateHetDelSegments,
        purity,
        ploidy);
    return TotalCopyNumberPrior::kDiploidBias;
  }

  // Build mBAF vector from the already-computed values
  vec<f64> candidate_mbafs_vec;
  candidate_mbafs_vec.reserve(candidate_segments.size());
  for (const auto& [seg, mbaf] : candidate_segments) {
    candidate_mbafs_vec.push_back(mbaf);
  }
  const arma::vec candidate_mbafs(candidate_mbafs_vec);

  // Fit KDE on per-segment mBAF values (already in [0, 0.5])
  const arma::vec weights(candidate_mbafs.n_elem, arma::fill::ones);
  const f64 wgd_bandwidth_multiplier = 1.1;
  const auto kde = GaussianKernelDensityEstimate(candidate_mbafs, weights, wgd_bandwidth_multiplier);

  auto local_maxima = FindLocalMaxima(kde.y);
  if (local_maxima.size() < kNumTopPeaks) {
    Logging::Info(
        "WGD QC: fewer than 2 peaks in clonal heterozygous deletion mBAF distribution, accepting purity/ploidy "
        "solution (purity={:.2f}, ploidy={:.2f})",
        purity,
        ploidy);
    return TotalCopyNumberPrior::kDiploidBias;
  }

  // Find the top 2 peaks by density
  std::partial_sort(local_maxima.begin(),
                    local_maxima.begin() + kNumTopPeaks,
                    local_maxima.end(),
                    [&kde](const size_t a, const size_t b) { return kde.y(a) > kde.y(b); });
  const f64 peak1 = kde.x(local_maxima[0]);
  const f64 peak2 = kde.x(local_maxima[1]);

  Logging::Info("WGD QC: top 2 clonal heterozygous deletion mBAF peaks at {:.3f} and {:.3f} (diff {:.3f})",
                peak1,
                peak2,
                std::abs(peak1 - peak2));

  if (std::abs(peak1 - peak2) > kWGDMaxPeakDifference) {
    Logging::Info(
        "WGD QC: peak difference exceeds {:.2f}, rejecting purity/ploidy solution (purity={:.2f}, ploidy={:.2f})",
        kWGDMaxPeakDifference,
        purity,
        ploidy);
    return TotalCopyNumberPrior::kUniform;
  }
  return TotalCopyNumberPrior::kDiploidBias;
}

/**
 * @brief Runs purity/ploidy search with an optional whole-genome-duplication (WGD) quality-control step.
 *
 * First performs a standard purity/ploidy search using a diploid-bias copy-number prior. If BAF data are available,
 * the result is then passed to EvaluatePurityPloidySolution. If that check determines the diploid solution is
 * inconsistent with the observed mBAF distribution (i.e. WGD is suspected), the search is re-run with a uniform
 * copy-number prior and the updated purity/ploidy estimates are returned instead.
 *
 * When BAF data are absent the WGD QC step is skipped and the diploid-bias result is returned directly.
 *
 * @param segments Genomic segments to use in the search.
 * @param logrs Per-locus log-ratio observations.
 * @param ref_obvs Per-locus reference allele depths (optional; required for WGD QC).
 * @param alt_obvs Per-locus alternate allele depths (optional; required for WGD QC).
 * @param options Purity/ploidy search options (seg_max_logr, seg_min_num_logrs, seg_min_num_snps).
 * @param grid_out_fname Path to write the purity/ploidy likelihood grid (diploid-bias run).
 * @param sex Sex of the sample (used to handle allosome segments correctly).
 * @param command_line_info Metadata to prepend to output files.
 * @return PurityPloidySearchOut containing the selected purity, ploidy, likelihood grids, and the
 *         TotalCopyNumberPrior that was ultimately used (kDiploidBias or kUniform).
 */
PurityPloidySearchOut PurityPloidySearchWithWGDQC(const std::vector<GenomicSegment>& segments,
                                                  const Observations& logrs,
                                                  const std::optional<Observations>& ref_obvs,
                                                  const std::optional<Observations>& alt_obvs,
                                                  const PurityPloidySearchOptions& options,
                                                  const fs::path& grid_out_fname,
                                                  const Sex sex,
                                                  const io::CommandLineInfo& command_line_info) {
  PurityPloidySearchOut out = PurityPloidySearch(segments,
                                                 logrs,
                                                 ref_obvs,
                                                 alt_obvs,
                                                 options.seg_max_logr,
                                                 options.seg_min_num_logrs,
                                                 options.seg_min_num_snps,
                                                 grid_out_fname,
                                                 sex,
                                                 command_line_info,
                                                 /*use_diploid_bias_prior=*/true);

  TotalCopyNumberPrior prior = TotalCopyNumberPrior::kDiploidBias;
  if (ref_obvs.has_value() && alt_obvs.has_value()) {
    prior = EvaluatePurityPloidySolution(
        segments, logrs, RefAltObservations{ref_obvs.value(), alt_obvs.value()}, out.purity, out.ploidy, options);
  }

  if (prior == TotalCopyNumberPrior::kUniform) {
    Logging::Info("WGD QC: re-running purity/ploidy search with uniform copy number prior");
    // Re-run with uniform prior, overwriting grid_out_fname so the persisted grid matches the chosen solution
    PurityPloidySearchOut uniform_out = PurityPloidySearch(segments,
                                                           logrs,
                                                           ref_obvs,
                                                           alt_obvs,
                                                           options.seg_max_logr,
                                                           options.seg_min_num_logrs,
                                                           options.seg_min_num_snps,
                                                           grid_out_fname,
                                                           sex,
                                                           command_line_info,
                                                           /*use_diploid_bias_prior=*/false);
    out.purity = uniform_out.purity;
    out.ploidy = uniform_out.ploidy;
    out.joint_ll_grid = std::move(uniform_out.joint_ll_grid);
    out.logr_ll_grid = std::move(uniform_out.logr_ll_grid);
    out.baf_ll_grid = std::move(uniform_out.baf_ll_grid);
    Logging::Info("WGD QC: uniform prior solution: purity={:.2f}, ploidy={:.2f}", out.purity, out.ploidy);
  }

  out.total_copy_number_prior = prior;
  return out;
}

}  // namespace xoos::cnc
