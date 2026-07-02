#include "likelihood/likelihood.h"

#include <cfloat>
#include <cmath>

#include <fmt/core.h>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>
#include <xoos/stats/copy-number-stats.h>
#include <xoos/stats/distributions.h>
#include <xoos/stats/stats.h>
#include <xoos/util/file-functions.h>
#include <xoos/util/math.h>

#include "copy-number-caller/copy-number-caller-options.h"
#include "io/copy-number-caller-default-filenames.h"
#include "io/write-segments.h"
#include "likelihood/baf-tools.h"
#include "likelihood/flag-copy-number-calls.h"
#include "mapq-utils.h"
#include "misc/sample-metadata-options.h"
#include "seg-to-vcf/seg-to-vcf.h"
#include "segmentation/read-segments.h"
#include "segmentation/segments-header.h"

namespace xoos::cnc {
using std::log;
using std::tuple;

const size_t kLikelihoodMaxCN2 = 20;  // for CN likelihood estimation after purity-ploidy search
constexpr size_t kNCN2 = kLikelihoodMaxCN2 + 1;

using segmentation::SegmentType;

// simple helper function to resolve sample ID for VCF output from the input segments
static std::string ResolveSampleIdFromSegments(const std::vector<GenomicSegment>& segments,
                                               const std::string& fallback_sample_id,
                                               const fs::path& segments_fname) {
  if (segments.empty() || fallback_sample_id != kSampleMetadataOptionsDefaultSampleId) {
    return fallback_sample_id;
  }
  for (const auto& seg : segments) {
    if (seg.id.has_value() && !seg.id->empty()) {
      return seg.id.value();
    }
  }
  return fallback_sample_id;
}

/**
 * @brief Use global tumor purity/ploidy estimates (user-provided or optionally estimate those here) to calculate the
 * joint logr-BAF likelihood of each segment. In calculating this joint likelihood we also return the most likely total
 * copy number of the segment as well as the the most likely minor copy number. See docs/likelihood.md for more
 * information
 * @param params
 */

void CalculateLikelihoodsMain(const CopyNumberCallerOptions& options) {
  using enum_util::FormatEnumName;
  const auto likelihood_out = (options.likelihood_options.mode == LikelihoodMode::kSomatic)
                                  ? options.output_dir / kDefaultSomaticCNCallsetSegOutput
                                  : options.output_dir / kDefaultGermlineCNCallsetSegOutput;
  const auto vcf_out = options.output_dir / kDefaultSomaticCNCallsetVcfOutput;
  std::vector<fs::path> output_paths = {likelihood_out};
  if (options.likelihood_options.mode == LikelihoodMode::kSomatic) {
    output_paths.push_back(vcf_out);
  }
  file::CheckFilePermissionsAndOutputPathExistence({}, output_paths);
  std::vector<GenomicSegment> segments =
      segmentation::ReadSegments(options.segments_fname.value(), options.likelihood_options.input_segments_type);
  const std::string sample_id =
      ResolveSampleIdFromSegments(segments, options.sample_metadata_options.sample_id, options.segments_fname.value());
  const segmentation::SegmentsHeader header = segmentation::ReadHeaderFromSegments(options.segments_fname.value());
  std::ifstream logr_ifstream(options.logrs_fname.value());
  Observations logrs = ReadObservations(logr_ifstream);
  LikelihoodModel ll_model = options.likelihood_options.likelihood_model;
  if ((ll_model == LikelihoodModel::kSerialSummarized) && !options.bafs_fname.has_value()) {
    Logging::Debug("likelihood model {} specified, but no BAFs provided. Defaulting to {}",
                   FormatEnumName(ll_model),
                   FormatEnumName(LikelihoodModel::kLogrSummarizedOnly));
    ll_model = LikelihoodModel::kLogrSummarizedOnly;
  }
  std::optional<Observations> ref_ads, alt_ads;
  if (options.bafs_fname.has_value()) {
    auto ref_alts = GetDepthsFromFile(options.bafs_fname.value());
    ref_ads = std::move(ref_alts.ref_obvs);
    alt_ads = std::move(ref_alts.alt_obvs);
  }
  std::optional<Observations> mapqs;
  if (options.mapping_qualities_fname.has_value()) {
    std::ifstream mapq_ifs(options.mapping_qualities_fname.value());
    mapqs = ReadObservations(mapq_ifs);
  }
  auto ll_opts = options.likelihood_options;
  segments = CalculateLikelihoods(
      segments, logrs, ref_ads, alt_ads, mapqs, ll_model, ll_opts, options.sample_metadata_options);
  for (auto& seg : segments) {
    seg.sex = options.sample_metadata_options.sex.value();
    seg.id = sample_id;
  }
  WriteSegments(likelihood_out, segments, options.likelihood_options.output_segments_type, options.command_line_info);
  if (options.likelihood_options.mode == LikelihoodMode::kSomatic) {
    const std::string reference_file =
        header.reference_file.empty() ? options.reference_genome_fname.string() : header.reference_file;
    const SeqLenMap seq_lengths = ParseFai(options.reference_genome_fai_fname);
    const std::unordered_map<std::string, size_t> seq_order = GetContigOrder(options.reference_genome_fai_fname);
    segmentation::WriteSegmentsToVcf(vcf_out,
                                     segments,
                                     {.seq_lengths = seq_lengths,
                                      .seq_order = seq_order,
                                      .sample_id = sample_id,
                                      .reference_file = reference_file,
                                      .sex = options.sample_metadata_options.sex.value(),
                                      .mode = CopyNumberCallerModes::kSomaticTumorNormalWGS,
                                      .purity = options.sample_metadata_options.purity.value(),
                                      .ploidy = options.sample_metadata_options.ploidy.value(),
                                      .vcf_purity_source = VcfPuritySource::kUserInput,
                                      .command_line_info = options.command_line_info});
  }
}

std::vector<GenomicSegment>& CalculateLikelihoods(std::vector<GenomicSegment>& segments,
                                                  const Observations& logrs,
                                                  const std::optional<Observations>& ref_ads,
                                                  const std::optional<Observations>& alt_ads,
                                                  const std::optional<Observations>& mapqs,
                                                  LikelihoodModel ll_model,
                                                  const LikelihoodOptions& options,
                                                  const SampleMetadataOptions& sample_metadata_options) {
  // calculate likelihoods
  if (!sample_metadata_options.sex.has_value()) {
    throw error::Error("Must provide sex to CalculateLikelihoods");
  }
  const Sex sex = sample_metadata_options.sex.value();
  if (!sample_metadata_options.purity.has_value() || !sample_metadata_options.ploidy.has_value()) {
    Logging::Error("Must provide purity and ploidy to CalculateLikelihoods!");
    throw std::runtime_error("purity and ploidy not provided to CalculateLikelihoods");
  }
  f64 purity = sample_metadata_options.purity.value();
  f64 ploidy = sample_metadata_options.ploidy.value();
  switch (ll_model) {
    case LikelihoodModel::kSerialSummarized: {
      if (!ref_ads.has_value() || !alt_ads.has_value()) {
        throw std::runtime_error("Serial model requires both ref and alt allele depths");
      }
      segments = AllSegmentsSerialLogRBAFLikelihoods(
          segments, logrs, ref_ads.value(), alt_ads.value(), purity, ploidy, options.mode, sex);
      break;
    }
    case LikelihoodModel::kLogrSummarizedOnly: {
      segments = AllSegmentsLogRLikelihoods(segments, logrs, purity, ploidy, options.mode, sex);
      break;
    }
    default: {
      Logging::Error("Invalid likelihood model specified for CalculateLikelihoods");
      throw std::runtime_error("no likelihood model specified (this should never happen)");
    }
  }
  // add metadata to the segments
  // also initiate an empty flags std::vector
  for (auto& seg : segments) {
    seg.sex = sample_metadata_options.sex.value();
    seg.id = sample_metadata_options.sample_id;
    seg.flags.emplace(std::vector<std::string>{});
  }
  // flag calls here
  // assign mean MAPQs and flag if necessary
  if (mapqs.has_value()) {
    arma::vec avg_mean_mapqs = GetAvgMeanMapqPerSegment(segments, mapqs.value());
    for (size_t i = 0; i < segments.size(); ++i) {
      segments[i].avg_mean_mapq = avg_mean_mapqs[i];
    }
    if (options.mapq_cutoff_for_calls > 0) {
      FlagCallsByMeanMapq(segments, options.mapq_cutoff_for_calls);
    }
  }
  // flag by expected LogR
  FlagCallsByExpectedTotalCopyNumber(segments, sex);
  FlagCallsByLength(segments, options.cnv_length_flag_min_size);
  FlagChrYCallsIfFemale(segments);
  return segments;
}

bool SkipBasedOnModeAndSex(const GenomicSegment& seg, LikelihoodMode mode, Sex sex) {
  // skip if:
  // 1. sex is unknown and in allosome
  // 2. sex is female and in chrY (keep all of chrX regardless of mode)
  // 3. sex is male and mode is somatic, and we are in chrY or a non-PAR region
  // 4. sex is male and mode is germline and in PAR region
  std::vector<bool> conditions{
      sex == Sex::kUnknown && seg.in_allosome.value_or(false),
      sex == Sex::kFemale && seg.in_allosome.value_or(false) && IsInChromY(seg.contig),
      sex == Sex::kMale && mode == LikelihoodMode::kSomatic && seg.in_allosome.value_or(false) &&
          (!seg.in_pseudo_autosomal_region.value_or(false) || IsInChromY(seg.contig)),
      sex == Sex::kMale && mode == LikelihoodMode::kGermline && seg.in_allosome.value_or(false) &&
          IsInChromY(seg.contig) && seg.in_pseudo_autosomal_region.value_or(false)};
  return std::any_of(conditions.begin(), conditions.end(), [](bool x) { return x; });
}

bool ExpectHaploid(const GenomicSegment& seg, LikelihoodMode mode, Sex sex) {
  // Only expect haploid in germline mode, in non-PAR region of an allosome, and if sex is male
  return (mode == LikelihoodMode::kGermline && seg.in_allosome.value_or(false) &&
          !seg.in_pseudo_autosomal_region.value_or(false) && sex == Sex::kMale);
}

/**
 * @brief Estimate total copy number per segment using only logr obserbations
 * @param segments
 * @param logrs
 * @param purity
 * @param ploidy
 * @return  std::vector of GenomicSegment
 */
std::vector<GenomicSegment>& AllSegmentsLogRLikelihoods(std::vector<GenomicSegment>& segments,
                                                        const Observations& logrs,
                                                        f64 purity,
                                                        f64 ploidy,
                                                        LikelihoodMode mode,
                                                        Sex sex) {
  IntervalTrees logr_trees(logrs);
  // we need to conisder all multiplicities from 0 to kMaxLikelihood2 inclusive
  arma::vec logr_likelihoods(kNCN2);  // idx=coverage
  f64 mean_segment_logr_sd = GetMeanSegmentLogRSD(segments, logrs, logr_trees);
  for (auto& seg : segments) {
    seg.ploidy = ploidy;
    seg.purity = purity;
    logr_likelihoods.fill(kNegativeInf);
    if (SkipBasedOnModeAndSex(seg, mode, sex)) {
      seg.total_copy_number = 0;
      seg.logr_likelihood = std::nullopt;
      seg.avg_mean_mapq = std::nullopt;
      seg.expected_total_copy_number = 0;
      seg.joint_likelihood = std::nullopt;
      continue;
    }
    bool expect_haploid = ExpectHaploid(seg, mode, sex);
    // If there are no observations and mean logR has no value, then this should represent a seed segment with no
    // overlapping observations. We can skip likelihood model predictions on this segment.
    if (seg.num_obs.value() == 0 && !seg.mean_logr.has_value()) {
      continue;
    }
    std::vector<size_t> logr_idxs = logr_trees.LookUp(seg.contig, seg.start, seg.end);
    Observations logr_obvs = logrs.FilterByIdxs(logr_idxs);
    SegmentAllSummarizedLogRLikelihoods(
        logr_obvs.obvs, mean_segment_logr_sd, purity, ploidy, expect_haploid, logr_likelihoods);
    // Choose best TCN state based on LogR likelihood alone
    size_t total_cn = logr_likelihoods.index_max();
    seg.total_copy_number = total_cn;
    seg.logr_likelihood = logr_likelihoods(total_cn);
    seg.joint_likelihood = logr_likelihoods(total_cn);
  }
  return segments;
}

/**
 * @brief Serial method for returning joint LogR-BAF likelihoods for all
 * segments. Here, we choose the total copy number state (TCN) solely by the
 * best LogR likelihood (as opposed to the joint method, which chooses the TCN
 * by the joint LogR-BAF likelihood).  The maximum BAF likelihood at this TCN
 * state is then added to the LogR likelihood at this state. The returned value
 * is this joint likelihood along with the corresponding TCN state and minor CN
 * state.
 * @param segments
 * @param logrs
 * @param ref_ads
 * @param alt_ads
 * @param purity
 * @param ploidy
 * @param mode (kGermline or kSomatic)
 * @param sex
 * @return  std::vector of GenomicSegment
 */
std::vector<GenomicSegment>& AllSegmentsSerialLogRBAFLikelihoods(std::vector<GenomicSegment>& segments,
                                                                 const Observations& logrs,
                                                                 const Observations& ref_ads,
                                                                 const Observations& alt_ads,
                                                                 f64 purity,
                                                                 f64 ploidy,
                                                                 LikelihoodMode mode,
                                                                 Sex sex) {
  IntervalTrees logr_trees(logrs);
  IntervalTrees variant_trees(ref_ads);
  f64 mean_segment_baf_sd = 0;
  mean_segment_baf_sd = GetMeanSegmentMBAFSD(segments, ref_ads, alt_ads, variant_trees);
  // we need to conisder all multiplicities from 0 to kMaxLikelihood2 inclusive
  arma::vec logr_likelihoods(kNCN2);  // idx=coverage
  arma::vec baf_likelihoods(kNCN2);   // row_idx=TCN, col_idx=multiplicity
  f64 mean_segment_logr_sd = GetMeanSegmentLogRSD(segments, logrs, logr_trees);
  for (auto& seg : segments) {
    seg.ploidy = ploidy;
    seg.purity = purity;
    if (SkipBasedOnModeAndSex(seg, mode, sex)) {
      seg.total_copy_number = 0;
      seg.minor_copy_number = 0;
      seg.major_copy_number = 0;
      seg.logr_likelihood = std::nullopt;
      seg.baf_likelihood = std::nullopt;
      seg.joint_likelihood = std::nullopt;
      seg.avg_mean_mapq = std::nullopt;
      seg.expected_total_copy_number = 0;
      continue;
    }
    bool expect_haploid = ExpectHaploid(seg, mode, sex);
    // If there are no observations and mean logR has no value, then this should represent a seed segment with no
    // overlapping observations. We can skip likelihood model predictions on this segment.
    if ((!seg.num_obs.has_value()) || (seg.num_obs.value() == 0 && !seg.mean_logr.has_value())) {
      continue;
    }
    // gather data from the segments
    std::vector<size_t> logr_idxs = logr_trees.LookUp(seg.contig, seg.start, seg.end);
    Observations logr_obvs = logrs.FilterByIdxs(logr_idxs);
    std::vector<size_t> variant_idxs = variant_trees.LookUp(seg.contig, seg.start, seg.end);
    Observations ref_ad_subset = ref_ads.FilterByIdxs(variant_idxs);
    Observations alt_ad_subset = alt_ads.FilterByIdxs(variant_idxs);
    if (logr_idxs.empty()) {
      Logging::Warn("Segment {}:{}-{} has 0 LogR observations. Probably a seed segment. Skipping...",
                    seg.contig,
                    seg.start + 1,
                    seg.end);
      continue;
    }
    if (variant_idxs.empty()) {
      Logging::Warn(
          "Segment {}:{}-{} has 0 BAF observations. Will set this segment's joint likelihood to its LogR likelihood, "
          "and set MajorCopyNumber and MinorCopyNumber to nan",
          seg.contig,
          seg.start + 1,
          seg.end);
      seg.num_snps = 0;
      seg.major_copy_number = std::nullopt;
      seg.minor_copy_number = std::nullopt;
      seg.mbaf = std::nullopt;
    } else {
      f64 mbaf = GetPeakOfBAFDistributionThenMirror(ref_ad_subset.obvs, alt_ad_subset.obvs);
      seg.mbaf = std::isnan(mbaf) ? std::nullopt : std::optional<f64>(mbaf);
    }
    // then, find all copy number likelihoods
    seg = SegmentSerialSummarizedLogRBAFLikelihood(seg,
                                                   logr_obvs.obvs,
                                                   purity,
                                                   ploidy,
                                                   mean_segment_logr_sd,
                                                   mean_segment_baf_sd,
                                                   false,
                                                   expect_haploid,
                                                   logr_likelihoods,
                                                   baf_likelihoods);
  }
  return segments;
}

const f64 kLikelihoodDefaultDF = 5;

f64 SegmentSummarizedLogRLikelihood(
    const arma::vec& logrs, f64 sd, f64 purity, f64 tumor_ploidy, f64 total_copy_number, bool expect_haploid) {
  f64 expected_logr = stats::ExpectedLogR(purity, tumor_ploidy, total_copy_number, expect_haploid);
  f64 mean_logr = arma::mean(logrs);
  f64 pdf = stats::NonStandardTPDF(mean_logr, kLikelihoodDefaultDF, expected_logr, sd);
  pdf = !math::IsCloseToZero(pdf) ? pdf : DBL_MIN;
  return log(pdf);
}

/**
 * @brief Calculates the BAF likelihood of a segment
 * @param mbaf - segment mirrored BAF estimate
 * @param sd - standard deviation of mirrored BAFs across all segments
 * @param purity - tumor purity
 * @param total_copy_number
 * @param multiplicity - minor copy number
 * @param expect_haploid
 * @return
 */
f64 SegmentSummarizedBAFLikelihood(
    f64 mbaf, f64 sd, f64 purity, f64 total_copy_number, f64 multiplicity, bool expect_haploid) {
  if (expect_haploid) {
    return NAN;
  }
  f64 expected_af = stats::ExpectedAF(purity, total_copy_number, multiplicity, expect_haploid);
  return log(stats::NonStandardTPDF(mbaf, kLikelihoodDefaultDF, expected_af, sd));
}

f64 GetMeanSegmentLogRSD(const std::vector<GenomicSegment>& segments,
                         const Observations& logrs,
                         IntervalTrees& logr_trees) {
  f64 total_sd = 0;
  f64 denom = 0;
  for (const auto& seg : segments) {
    if (IsInAllosome(seg.contig)) {
      continue;
    }
    std::vector<size_t> logr_idxs = logr_trees.LookUp(seg.contig, seg.start, seg.end);
    Observations logr_obvs = logrs.FilterByIdxs(logr_idxs);
    f64 sd = arma::stddev(logr_obvs.obvs);
    if (std::isnan(sd)) {
      throw error::Error("Segment {}:{}-{} has NaN logR standard deviation", seg.contig, seg.start, seg.end);
    }
    total_sd += sd;
    denom += 1;
  }

  if (math::IsCloseToZero(denom)) {
    throw error::Error("No segments on autosomal chromosomes; cannot calculate mean segment logR standard deviation");
  }

  return total_sd / denom;
}

/**
 * @brief get the standard deviation of the weighted mBAF means of a set of segments
 * @param segments
 * @param ref_ads
 * @param alt_ads
 * @param variant_trees
 * @return
 */
f64 GetMeanSegmentMBAFSD(const std::vector<GenomicSegment>& segments,
                         const Observations& ref_ads,
                         const Observations& alt_ads,
                         IntervalTrees& variant_trees) {
  f64 num = 0;
  f64 denom = 0;
  for (const auto& seg : segments) {
    if (IsInAllosome(seg.contig)) {
      continue;
    }
    std::vector<size_t> variant_idxs = variant_trees.LookUp(seg.contig, seg.start, seg.end);
    if (variant_idxs.size() > 1) {
      Observations ref_ad_subset = ref_ads.FilterByIdxs(variant_idxs);
      Observations alt_ad_subset = alt_ads.FilterByIdxs(variant_idxs);
      f64 mbaf_mean = GetPeakOfBAFDistributionThenMirror(ref_ad_subset.obvs, alt_ad_subset.obvs);
      num += WeightedMBAFSD(ref_ad_subset.obvs, alt_ad_subset.obvs, mbaf_mean);
      denom += 1;
    }
  }

  if (math::IsCloseToZero(denom)) {
    throw error::Error("No autosomal segments with >1 variant; cannot calculate mean segment mBAF standard deviation");
  }

  return num / denom;
}

/**
 * @brief Get the summarized likelihoods (mean of point for all total CNs from 0 -> ret.size()
 * likelihood for TCN=c: likelihood(mean(logrs) | c)
 * @param logrs - std::vector of observations associated with the segment
 * @param sd - standard deviation of mean logr across all segments
 * @param purity - tumor purity
 * @param tumor_ploidy
 * @param ret - std::vector to fill. size of std::vector determines how many likelihoods to calculate [0, ret.size()-1]
 * @return reference to ret
 */
arma::vec& SegmentAllSummarizedLogRLikelihoods(
    const arma::vec& logrs, f64 sd, f64 purity, f64 tumor_ploidy, bool expect_haploid, arma::vec& ret) {
  size_t c = 0;
  ret.for_each([&logrs, sd, purity, tumor_ploidy, &c, expect_haploid](arma::vec::elem_type& x) {
    x = SegmentSummarizedLogRLikelihood(logrs, sd, purity, tumor_ploidy, static_cast<f64>(c++), expect_haploid);
  });
  return ret;
}

/**
 * @brief Get all the summarized BAF likelihoods associated with a segment. The number of likelihoods is each possible
 * minor copy number, ie. total_copy_number / 2 likelihood for TCN=c, MinorCN=m: likelihood(mean(baf) | c,m) for m in
 * [min(a,d) for a,d in zip(ref_ads, alt_ads)]
 * @param ref_ads - reference AD observations associated with the segment
 * @param alt_ads - alt AD observations associated with the segment
 * @param sd - standard deviation of mean mBAF across all segments
 * @param purity  - tumor purity
 * @param total_copy_number
 * @param ret - std::vector to fill. size of std::vector must be > total_copy_number/2
 * @return reference to ret
 */
SegmentAllSummarizedBAFLikelihoodsResults SegmentAllSummarizedBAFLikelihoods(
    f64 mbaf, f64 sd, f64 purity, f64 total_copy_number, bool expect_haploid) {
  // calculate baf likelihood for each TCN/MinorCN combo
  auto tcn_int = static_cast<size_t>(total_copy_number);
  arma::vec baf_likelihoods(kNCN2);  // row_idx=TCN, col_idx=multiplicity
  baf_likelihoods.fill(kNegativeInf);
  for (size_t m = 0; m <= tcn_int / 2; ++m) {
    baf_likelihoods(m) =
        SegmentSummarizedBAFLikelihood(mbaf, sd, purity, total_copy_number, static_cast<f64>(m), expect_haploid);
  }
  return {baf_likelihoods};
}

static f64 GetTcnPrior(size_t n_cns, size_t cn) {
  // Consider the copy neutral state (2) as 2 possible states. This makes the total number of states n_cns + 1.
  // This allows the copy neutral state to have f64 the prior.
  auto n = static_cast<f64>(n_cns) + 1;
  return cn == 2 ? log(2.0 / n) : log(1.0 / n);
}

GenomicSegment& SegmentSerialSummarizedLogRBAFLikelihood(GenomicSegment& segment,
                                                         const arma::vec& logrs,
                                                         f64 purity,
                                                         f64 ploidy,
                                                         f64 mean_logr_sd,
                                                         f64 mean_baf_sd,
                                                         bool prior_2n,
                                                         bool expect_haploid,
                                                         arma::vec& logr_likelihoods,
                                                         arma::vec& baf_likelihoods) {
  logr_likelihoods.fill(kNegativeInf);
  SegmentAllSummarizedLogRLikelihoods(logrs, mean_logr_sd, purity, ploidy, expect_haploid, logr_likelihoods);
  // Choose best TCN state based on LogR likelihood alone
  if (prior_2n) {
    for (size_t i = 0; i < logr_likelihoods.size(); ++i) {
      logr_likelihoods[i] += GetTcnPrior(logr_likelihoods.size(), i);
    }
  }
  size_t total_cn = logr_likelihoods.index_max();
  f64 logr_likelihood = logr_likelihoods(total_cn);
  // calculate baf likelihoods for the total CN determined from the logr_likelihoods, over all minor CNs
  if (segment.mbaf.has_value() && !expect_haploid) {
    auto baf_results = SegmentAllSummarizedBAFLikelihoods(
        segment.mbaf.value(), mean_baf_sd, purity, static_cast<f64>(total_cn), expect_haploid);
    baf_likelihoods = baf_results.baf_likelihoods;
    size_t minor_cn = baf_likelihoods.index_max();
    segment.total_copy_number = total_cn;
    segment.minor_copy_number = minor_cn;
    segment.major_copy_number = total_cn - minor_cn;
    segment.logr_likelihood = logr_likelihood;
    segment.baf_likelihood = baf_likelihoods(minor_cn);
    segment.joint_likelihood = segment.logr_likelihood.value() + segment.baf_likelihood.value();

  } else {
    segment.total_copy_number = total_cn;
    segment.logr_likelihood = logr_likelihood;
    segment.baf_likelihood = std::nullopt;
    segment.joint_likelihood = logr_likelihood;
  }
  segment.ploidy = ploidy;
  segment.purity = purity;
  return segment;
}

}  // namespace xoos::cnc
