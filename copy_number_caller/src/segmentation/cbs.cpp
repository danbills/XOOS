/**
 * reference:
 * Olshen et al, 2004. "Circular binary segmentation for the analysis of array‐based DNA copy number data"
 * https://academic.oup.com/biostatistics/article/5/4/557/275197
 */

#include "segmentation/cbs.h"

#include <cmath>
#include <limits>
#include <ranges>  // NOLINT cpplint doesn't yet support some C++20 headers
#include <stack>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/stats/external/asa243/asa243.hpp>
#include <xoos/util/math.h>

#include "segmentation/boundaries.h"
#include "segmentation/cbs-util.h"
#include "segmentation/max-t.h"

namespace xoos::cnc::segmentation {

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::chrono::nanoseconds;
using std::chrono::seconds;

const size_t kCbsSiegmundPIters = 100;
const size_t kCbsMinObvsForRecursion = 10;
const size_t kCbsHybridK = 25;
const size_t kCbsMinObvsForShortcut = 10;
const f64 kCbsMaxTForAutomaticNonRejectNull = 0.1;
const size_t kCbsMaxStackDepth = 1024;
const size_t kCbsDefaultNPermutations = 10000;
const f64 kCbsDefaultMaxP = 0.01;

/**
 * @brief permute null distribution describing scenario where the maximal
 * t-value of a circular parition of the observation std::vector is not significantly
 * deviant from the maximal t-value of a circular partitioning of any random
 * permutation of the observations
 * acccept null if (max_p*n_permutations) permutations' max-t statistics exceed t
 * reject null if above condition fails, or if boundaries are provided and the boundary test fails at any permutation
 * @param obvs observation std::vector of doubles
 * @param n_permutations max number of permutations to test
 * @param max_p maximum p-value for rejecting null
 * @param t maximal t-statistic for circular partitioning of observation std::vector
 * @param max_obvs_size maximum number of observations allowed per segment
 * @param boundaries pre-calculated boundaries for early-stopping
 * @param min_obs_per_segment minimum number of observations allowed per segment
 * @return true if null is rejected else false
 */
bool CircularPermuteNullAndTestSignificance(const arma::vec& obvs,
                                            size_t n_permutations,
                                            f64 max_p,
                                            Segment max_t_seg,
                                            size_t max_obs_size,
                                            const arma::uvec& boundaries,
                                            size_t min_obs_per_segment,
                                            CbsMaxTMethod max_t_method) {
  auto max_n_gt_t = static_cast<size_t>(max_p * static_cast<f64>(n_permutations));
  Logging::Debug("Estimating null over maximum {} iters using early stopping measure. hard stop at {} excds",
                 n_permutations,
                 max_n_gt_t + 1);
  size_t n_gt_t = 0;
  // k denotes start-index in boundaries std::vector
  size_t k = max_n_gt_t * (max_n_gt_t + 1) / 2;
  arma::vec shuffled_obvs(obvs);
  Logging::Debug("boundary (std::vector size {}) start index is: {}", boundaries.n_elem, k);
  arma::vec centered_obvs = obvs - arma::mean(obvs);
  f64 tss = arma::sum<arma::vec>(arma::pow(centered_obvs, 2));
  for (size_t i = 0; i < n_permutations; ++i) {
    arma::arma_rng::set_seed(i);
    FisherYatesShuffleInPlace(shuffled_obvs);
    Segment perm_max_t_seg = GetMaxTBreakpoints(shuffled_obvs, max_obs_size, min_obs_per_segment, max_t_method);
    f64 new_tss = tss < perm_max_t_seg.t ? perm_max_t_seg.t + 1 : tss;
    perm_max_t_seg.t = perm_max_t_seg.t / ((new_tss - perm_max_t_seg.t) / (static_cast<f64>(obvs.size()) - 2));
    if (max_t_seg.t < perm_max_t_seg.t) {
      n_gt_t += 1;
    }
    // if the number of non-reject permutations exceeds <max_p*n_permutations>, then we can automatically say that the
    // segment isn't significant
    if (n_gt_t > max_n_gt_t) {
      Logging::Debug(
          "stopped permuting for segment of size {} at {} iters {} excds. no reject null", obvs.n_elem, i, n_gt_t);
      return false;
    }
    if (boundaries.n_elem > 0 && boundaries[k + n_gt_t] < i) {
      Logging::Debug("stopped permuting for segment of size {} at {} iters. reject null", obvs.n_elem, i);
      return true;
    }
  }
  Logging::Debug("permuted maximum amount of times for segment of size {} at {} iters", obvs.n_elem, n_permutations);
  return static_cast<f64>(n_gt_t) / static_cast<f64>(n_permutations) < max_p;
}

/**
 * @brief permute null distribution describing scenario where the t-value of a
 * binary parition of the observation std::vector is not significantly deviant from
 * the maximal t-value of a proportionally similar binary partitioning of any
 * random permutation of the observations
 * acccept null if (max_p*n_permutations) permutations' max-t statistics exceed t
 * reject null if above condition fails
 * @param obvs observation std::vector of doubles
 * @param n_permutations max number of permutations to test
 * @param max_p maximum p-value for rejecting null
 * @param seg1 first segment in the binary partition
 * @param seg2 second segment in the binary partition
 * @return true if null is rejected else false
 * TODO: pass in size_ts instead of Segments
 */
bool BinaryPermuteNullAndTestSignificance(const arma::vec& obvs,
                                          size_t n_permutations,
                                          f64 max_p,
                                          const Segment& seg1,
                                          const Segment& seg2,
                                          f64 min_t_for_automatic_reject_null) {
  arma::vec v1 = obvs(arma::span(seg1.start, seg1.end - 1));
  arma::vec v2 = obvs(arma::span(seg2.start, seg2.end - 1));
  f64 t = std::abs(TStatFromVecs(v1, v2));
  // shortcut: reject null if both segments are large enough and t is large
  if (t > min_t_for_automatic_reject_null && v1.size() > kCbsMinObvsForShortcut && v2.size() > kCbsMinObvsForShortcut) {
    return true;
  }
  // otherwise if both segments are large enough, we can use the student's t distribution to assess if both segments are
  // signiciantly different from each other
  if (v1.size() > kCbsHybridK && v2.size() > kCbsHybridK) {
    auto df = static_cast<f64>(v1.size() + v2.size() - 2);
    f64 p = 2 * stats::tnc(-std::abs(t), df, 0);
    return p < max_p;
    // if one or both segments are very small, then we have to use permutation testing because it's not enough
    // observations to assume that the data follows a normal distribution
  } else {
    arma::vec to_cmp = arma::join_cols(v1, v2);
    f64 total_sum = arma::sum(to_cmp);
    f64 total_sum_of_squares = arma::sum<arma::vec>(arma::pow(to_cmp, 2));
    auto minor_arc_len = static_cast<f64>(std::min(v1.size(), v2.size()));
    auto major_arc_len = static_cast<f64>(std::max(v1.size(), v2.size()));
    f64 nexcd = 0.0;
    f64 max_excd = static_cast<f64>(n_permutations) * max_p + 1;
    for (size_t i = 0; i < n_permutations; ++i) {
      // TODO make sure n_permutations is less than the total number of possible permutations!!
      // since we know the minor arc is small, we don't need to permute the
      // entire array. We can just sample minor_arc_len points, calculate sums for those directly and infer the sums for
      // the major arc from the minor arc sum
      arma::arma_rng::set_seed(i);
      auto minor_arc_sample_idxs(arma::randi<arma::uvec>(static_cast<size_t>(minor_arc_len),
                                                         arma::distr_param(0, static_cast<s32>(to_cmp.size() - 1))));
      f64 minor_arc_sum = 0;
      f64 minor_arc_sum_of_squares = 0;
      for (size_t idx : minor_arc_sample_idxs) {
        minor_arc_sum += to_cmp(idx);
        minor_arc_sum_of_squares += std::pow(to_cmp(idx), 2);
      }
      f64 major_arc_sum = total_sum - minor_arc_sum;
      f64 major_arc_sum_of_squares = total_sum_of_squares - minor_arc_sum_of_squares;
      // calculate the t statistic
      f64 minor_arc_mean = minor_arc_sum / minor_arc_len;
      f64 minor_arc_var = Variance(minor_arc_len, minor_arc_mean, minor_arc_sum_of_squares);
      f64 major_arc_mean = major_arc_sum / major_arc_len;
      f64 major_arc_var = Variance(major_arc_len, major_arc_mean, major_arc_sum_of_squares);
      f64 perm_t =
          std::abs(TStat(minor_arc_len, minor_arc_mean, minor_arc_var, major_arc_len, major_arc_mean, major_arc_var));
      if (perm_t > t) {
        nexcd += 1;
      }
      if (nexcd >= max_excd) {
        Logging::Debug("adjacent segs of size {} / mean {} and size {} / mean {} are not significantly different, p={}",
                       seg1.end - seg1.start,
                       arma::mean(obvs(arma::span(seg1.start, seg1.end - 1))),
                       seg2.end - seg2.start,
                       arma::mean(obvs(arma::span(seg2.start, seg2.end - 1))),
                       static_cast<f64>(nexcd) / static_cast<f64>(n_permutations));
        return false;
      }
    }
    return static_cast<f64>(nexcd) / static_cast<f64>(n_permutations) < max_p;
  }
}

const f64 kCbsSiegmundPMinChange = 1e6;
const f64 kCbsSiegmundPMaxIter = 100000;

/**
 * @brief described in Siegmund, "Approximate Tail Probabilities for the Maxima of Some Random Fields", 1988
 * @param x the parameter for Siegmund apprxomiation Nu function
 * @return the "nu" parameter for the Siegment p-value approximation equation
 *
 */
f64 SiegmundPNuFull(f64 x) {
  auto fn = [&x](f64 l) {
    f64 power = std::pow(l, -1);
    f64 half_x = -0.5 * x;
    f64 sqrt_l = std::pow(l, 0.5);
    arma::vec vec = arma::vec{half_x * sqrt_l};
    return power * arma::normcdf(vec)[0];
  };
  auto change_frac = [](const f64 p, const f64 n) {
    return math::IsCloseToZero(p) ? std::numeric_limits<f64>::max() : std::abs(p - n) / p;
  };
  // need to initialize ptotal and total to avoid div by 0
  f64 ptotal = fn(1);
  f64 total = ptotal + fn(2);
  f64 l = 3;
  while (change_frac(ptotal, total) > kCbsSiegmundPMinChange && l < kCbsSiegmundPMaxIter) {
    ptotal = total;
    total += fn(l);
    l += 1;
  }
  f64 ret = 2 * std::pow(x, -2);
  ret *= std::exp(-2 * total);
  return ret;
}

/**
 * @brief P-value approximation described in Siegmund, "Approximate Tail
 * Probabilities for the Maxima of Some Random Fields", 1986. This is used for
 * when you want a p-value from comparing the maxima of two samples. For
 * intuition, if you're comparing the means instead of the maxima then you can
 * get away with just using the normal CDF.
 * @param t test statistic for which to calculate significance
 * @param m size of sample
 * @param k scale parameter
 */
f64 SiegmundP(f64 t, f64 m, f64 k) {
  f64 pdf_t = arma::normpdf(arma::vec{t})[0];
  f64 p = 0.25 * std::pow(t, 3) * pdf_t;
  f64 delta = k / m;
  f64 incr = (1 - delta - 0.5) / static_cast<f64>(kCbsSiegmundPIters);
  arma::vec fn_to_intg_x(kCbsSiegmundPIters);
  arma::vec fn_to_intg_y(kCbsSiegmundPIters);
  for (size_t i = 0; i < kCbsSiegmundPIters; ++i) {
    f64 k = 0.5 + (static_cast<f64>(i) * incr);
    f64 k_1_min_k = k * (1 - k);
    f64 nu_sq = std::pow(SiegmundPNuFull(t / std::pow(m * k_1_min_k, 0.5)), 2);
    f64 y = nu_sq / std::pow(k_1_min_k, 2);
    fn_to_intg_x[i] = k;
    fn_to_intg_y[i] = y;
  }
  arma::mat integral = arma::trapz(fn_to_intg_x, fn_to_intg_y);
  if (integral.n_elem != 1) {
    throw error::Error("cbs: SiegmundP trapz returned unexpected number of elements");
  }
  f64 integral_area = integral[0];
  p = integral_area * p;
  return p;
}

/**
 * @brief the "hybrid" test for significance of a circular segmentation. If the
 * two segments are large enough, then use Siegmund's approximation to calculate
 * the portion of the p-value associated with >k sized segments. If the
 * partial-p is significant, use permutation testing to calculate the portion of
 * the p-value associated with <k sized segments and make a final decision on
 * significance. If the one of the segments is too small, default to permutation
 * testing for the whole p-value.
 * @param obvs observation std::vector
 * @param n_permutations max number of permutations to test
 * @param max_p maximum p-value for rejecting null
 * @param t maximal t-statistic for circular partitioning of observation std::vector
 * @param boundaries pre-calculated boundaries for early-stopping
 * @param min_obs_per_segment minimum number of observations allowed per segment
 * @return true if the segmentation was significant according to the p-value calculation, false if not
 */
bool HybridTestForSignificance(const arma::vec& obvs,
                               size_t n_permutations,
                               f64 max_p,
                               Segment max_t_seg,
                               const arma::uvec& boundaries,
                               size_t min_obs_per_segment,
                               CbsMaxTMethod max_t_method,
                               f64 min_t_for_automatic_reject_null) {
  // if the two arcs are large enough, then we can reject null if t is a large-ish number
  if (max_t_seg.end - max_t_seg.start > kCbsMinObvsForShortcut &&
      obvs.size() - max_t_seg.end + max_t_seg.start > kCbsMinObvsForShortcut &&
      max_t_seg.t > min_t_for_automatic_reject_null) {
    Logging::Debug("{} elem segment because t ({}) > {}", obvs.n_elem, max_t_seg.t, min_t_for_automatic_reject_null);
    return true;
    // we can also not-reject null automaticaly if t is realy small
  } else if (max_t_seg.t < kCbsMaxTForAutomaticNonRejectNull) {
    Logging::Debug("not segmenting {} elem segment because t ({}) < {}",
                   obvs.n_elem,
                   max_t_seg.t,
                   kCbsMaxTForAutomaticNonRejectNull);
    return false;
  }
  // if there aren't enough observations, we cannnot calculate approximate p. We have to do permutation testing on the
  // whole thing
  if (obvs.n_elem < kCbsHybridK * 4) {
    Logging::Debug("{} elem segment: using default P method", obvs.n_elem);
    return CircularPermuteNullAndTestSignificance(
        obvs, n_permutations, max_p, max_t_seg, 0, boundaries, min_obs_per_segment, max_t_method);
  } else {
    // if there enough observations, we can approximate p for all the partitions where both arcs are big enough
    Logging::Debug("{} elem segment: using approx P method", obvs.n_elem);
    // get null distribution for all splits minor_arc>=k && major_arc>=k
    f64 p = SiegmundP(max_t_seg.t, static_cast<f64>(obvs.n_elem), static_cast<f64>(kCbsHybridK));
    // if this is significant, then we calculate the part of the p-value for all the paritions where one arc is small
    if (p < max_p) {
      Logging::Debug("{} elem segment: approx p is significant ({}), will also run P method", obvs.n_elem, p);
      if (max_t_method == CbsMaxTMethod::kFast) {
        return CircularPermuteNullAndTestSignificance(obvs,
                                                      n_permutations,
                                                      max_p - p,
                                                      max_t_seg,
                                                      kCbsHybridK,
                                                      boundaries,
                                                      min_obs_per_segment,
                                                      CbsMaxTMethod::kFastMaxN);
      } else {
        return CircularPermuteNullAndTestSignificance(obvs,
                                                      n_permutations,
                                                      max_p - p,
                                                      max_t_seg,
                                                      kCbsHybridK,
                                                      boundaries,
                                                      min_obs_per_segment,
                                                      CbsMaxTMethod::kBruteForce);
      }
    }
  }
  return false;
}

/**
 * @brief Find the best scoring sub-segment (by t-statistic) and, if
 * statistically significant, return all segments created from this best-scoring
 * sub-segment's breakpoints. If not significant, return an empty std::vector
 * @param obvs observation std::vector
 * @param max_p maximum p-value for rejecting null
 * @param n_permutations max number of permutations to test
 * @param boundaries pre-calculated boundaries for early-stopping
 * @param min_obs_per_segment minimum number of observations allowed per segment
 * @return a maximum of three segments w.r.t to obvs. These segments should span the entirety of obvs
 */
std::vector<Segment> GetSegments(const arma::vec& obvs,
                                 f64 max_p,
                                 size_t n_permutations,
                                 const arma::uvec& boundaries,
                                 size_t min_obs_per_segment,
                                 CbsMaxTMethod max_t_method,
                                 f64 min_t_for_automatic_reject_null) {
  // get the highest t-stat and corresponding segment for the real observations
  Logging::Debug("finding max-T sub-segment for {} obvs segment", obvs.n_elem);
  arma::vec centered_obvs = obvs - arma::mean(obvs);
  f64 tss = arma::sum<arma::vec>(arma::pow(centered_obvs, 2));
  Segment max_t_seg = GetMaxTBreakpoints(centered_obvs, 0, min_obs_per_segment, max_t_method);
  tss = tss < max_t_seg.t ? max_t_seg.t + 1 : tss;
  max_t_seg.t = max_t_seg.t / ((tss - max_t_seg.t) / (static_cast<f64>(obvs.size()) - 2));
  Logging::Debug("found max-T sub-segment for {} obvs segment", obvs.n_elem);
  // conditions for immediate stop on recursion (no change, segmentation too small)
  if (math::IsCloseToZero(max_t_seg.t)) {
    Logging::Debug("segment with {} obvs is completely uniform", obvs.n_elem);
    return {};
  }
  if (max_t_seg.end - max_t_seg.start < min_obs_per_segment) {
    Logging::Debug("segment with {} obvs has very small max-t sub-segment", obvs.n_elem);
    return {};
  }
  if (!HybridTestForSignificance(obvs,
                                 n_permutations,
                                 max_p,
                                 max_t_seg,
                                 boundaries,
                                 min_obs_per_segment,
                                 max_t_method,
                                 min_t_for_automatic_reject_null)) {
    return std::vector<Segment>{};
  }
  Logging::Debug("found significance for max-T sub-segment for {} obvs segment", obvs.n_elem);
  std::vector<Segment> out_segments;
  if (max_t_seg.start != 0) {
    out_segments.emplace_back(0, max_t_seg.start);
  }
  out_segments.emplace_back(max_t_seg.start, max_t_seg.end);
  if (max_t_seg.end != obvs.n_elem) {
    out_segments.emplace_back(max_t_seg.end, obvs.n_elem);
  }
  return out_segments;
}

arma::vec TruncateOutliers(const arma::vec& obvs) {
  arma::vec ret(obvs);
  // get the 0.1 and .99
  arma::vec quantiles = arma::quantile(obvs, arma::vec({0.01, 0.99}));
  f64 interquantile_range = quantiles[1] - quantiles[0];
  f64 lower_bound = quantiles[0] - 1.5 * interquantile_range;
  f64 upper_bound = quantiles[1] + 1.5 * interquantile_range;
  for (auto& x : ret) {
    // truncate if the observation is above or below the quantiles
    if (x < lower_bound) {
      x = lower_bound;
    }
    if (x > upper_bound) {
      x = upper_bound;
    }
  }
  return ret;
}

/**
 * @brief Segment an observation std::vector using the Circular Binary Segmentation
 * Algorithm described in Olshen et al., 2004, "Circular binary segmentation
 * for the analysis of array- based DNA copy number data", 2004.
 * and implementing the speed-ups described in Venkatraman and Olshen, 2007,
 * "A faster circular binary segmentation algorithm for the analysis of array
 * CGH data", 2007.
 * @param obvs observation std::vector
 * @param max_p maximum p-value for rejecting null
 * @param n_permutations max number of permutations to test
 * @param min_obs_per_segment minimum number of observations allowed per segment
 * @return all significant segments found in obvs
 */
std::vector<Segment> CircularBinarySegmentation(const arma::vec& obvs, const SegmentationOptions& options) {
  size_t n = obvs.n_elem;
  // we need to truncate outliers to the same max value because one or two very extreme outliers can cause significance
  // testing to fail, causing severe under-segmentation
  arma::vec use_these_obvs(obvs);
  if (options.truncate_outliers) {
    Logging::Info("truncating outliers");
    use_these_obvs = TruncateOutliers(obvs);
  }
  auto time1 = high_resolution_clock::now();
  if (n == 0) {
    Logging::Info("no observations to segment!");
    return {};
  }
  std::vector<Segment> ret;
  std::stack<Segment> seg_stack;
  Segment first_seg = {.start = 0, .end = n};
  seg_stack.push(first_seg);
  while (!seg_stack.empty()) {
    Logging::Debug("stack has {} elements", seg_stack.size());
    const Segment& segment = seg_stack.top();
    seg_stack.pop();
    // If segment is small enough that we can't perform permutation testing if we break it up, then don't recurse
    if (segment.end - segment.start < kCbsMinObvsForRecursion) {
      Logging::Debug("adding small segment with {} obvs", (segment.end - segment.start));
      ret.push_back(segment);
    } else {
      // NOTE: we can probably make this faster if we pass the indices into GetSegments instead of copying
      arma::vec sub_obvs = use_these_obvs(arma::span(segment.start, segment.end - 1));
      // "recursive" step - partition segment into sub-segments
      arma::uvec boundaries;
      if (options.n_permutations == kCbsDefaultNPermutations &&
          math::IsEqualWithTolerance(options.max_p, kCbsDefaultMaxP)) {
        boundaries = kBoundaryP10000A01Eta05;
      }
      Logging::Debug("Segmenting on {} obs", sub_obvs.n_elem);
      std::vector<Segment> sub_segments = GetSegments(sub_obvs,
                                                      options.max_p,
                                                      options.n_permutations,
                                                      boundaries,
                                                      options.min_obs_per_segment,
                                                      options.cbs_method,
                                                      options.min_t_for_automatic_segmentation);
      Logging::Debug("Segmented on {} obs", sub_obvs.n_elem);
      // merge back single-obs sub-segments
      if (options.no_single_obvs_subsegments) {
        sub_segments = MergeSingleObvsEndSegment(sub_segments, options.min_obs_per_segment);
      }
      if (sub_segments.size() <= 1) {
        // if this segment cannot be broken down, then return it as one segment
        Logging::Debug("adding un-breakable segment with {} obvs", (segment.end - segment.start));
        ret.push_back(segment);
      } else if (!options.disable_merging) {
        Logging::Debug("considering whether to merge back segments");
        // "verify" that the segments are actually different from each other by doing the binary segmentation t-test
        // procedure .
        // First push all segments onto a stack in reverse order, such that the first segment is on the top of the
        // stack, etc.
        std::stack<Segment> segment_merge_stack;
        for (auto it : sub_segments | std::views::reverse) {
          segment_merge_stack.push(it);
        }
        sub_segments.clear();
        // while there are at least 2 elements in the stack, compare the top two
        // Segments in the stack. If they fail the binary segmentation test of
        // significance, then merge these two segments, and add the merged
        // segment onto the stack. If they are significantly different, then add
        // the top Segment to the array, and keep the second Segment in the
        // stack, so that it can be tested against the next Segment in the
        // following iteration
        // The loop necessarily stops when there is only 1 item in the stack.
        // This item is then added to the return std::vector
        while (segment_merge_stack.size() > 1) {
          Segment seg_to_merge1 = segment_merge_stack.top();
          segment_merge_stack.pop();
          Segment seg_to_merge2 = segment_merge_stack.top();
          if (!BinaryPermuteNullAndTestSignificance(sub_obvs,
                                                    options.n_permutations,
                                                    options.max_p,
                                                    seg_to_merge1,
                                                    seg_to_merge2,
                                                    options.min_t_for_automatic_segmentation)) {
            Segment new_sub_segment(seg_to_merge1.start, seg_to_merge2.end);
            segment_merge_stack.pop();  // remove seg_to_merge2 from top
            segment_merge_stack.push(new_sub_segment);
          } else {
            sub_segments.push_back(seg_to_merge1);
            // seg_to_merge2 is kept at the top of the stack
          }
        }
        sub_segments.push_back(segment_merge_stack.top());
        // case for when *all* the subsegments ended up being merged back together
        if (sub_segments.size() == 1) {
          Logging::Debug("adding un-breakable segment with {} obvs", (segment.end - segment.start));
          ret.push_back(segment);
        } else {
          // these sub-segments make no assumption about which "parent" its observations come from. Therefore, we have
          // to re-adjust the start/end fields of the sub-segments with respect to the parent std::vector's coordinate
          // system
          for (auto new_segment : sub_segments) {
            new_segment.start += segment.start;
            new_segment.end += segment.start;
            if (new_segment.end > n) {
              throw error::Error("cbs: sub-segment end exceeds observation count after coordinate adjustment");
            }
            Logging::Debug("pushing segment with {} obvs", (new_segment.end - new_segment.start));
            seg_stack.push(new_segment);
            if (seg_stack.size() > kCbsMaxStackDepth) {
              throw std::runtime_error("cbs: max depth exceeded!");
            }
          }
        }
      } else {
        for (auto new_segment : sub_segments) {
          new_segment.start += segment.start;
          new_segment.end += segment.start;
          if (new_segment.end > n) {
            throw error::Error("cbs: sub-segment end exceeds observation count after coordinate adjustment");
          }
          Logging::Debug("pushing segment with {} obvs", (new_segment.end - new_segment.start));
          seg_stack.push(new_segment);
          if (seg_stack.size() > kCbsMaxStackDepth) {
            throw std::runtime_error("cbs: max depth exceeded!");
          }
        }
      }
    }
  }
  std::sort(ret.begin(), ret.end(), [](const Segment& l, const Segment& r) { return l.start < r.start; });
  auto time2 = high_resolution_clock::now();
  auto total_time = duration_cast<milliseconds>(time2 - time1).count();
  Logging::Debug("found {} segments over {} elements in {}ms", ret.size(), n, total_time);
  return ret;
}

/**
 * @brief Given a set of 2-3 subsegments, if the first or the last sub-segment has fewer than min_obs_per_segment
 * observations, then merge it with its adjacent sub-segment
 * @param sub_segments a std::vector of 2-3 sub-segments
 * @param min_obs_per_segment  minimum number of observations allowed per segment
 * @return
 */
std::vector<Segment> MergeSingleObvsEndSegment(const std::vector<Segment>& sub_segments, size_t min_obs_per_segment) {
  if (sub_segments.size() <= 1) {
    return sub_segments;
  }
  std::vector<Segment> res;
  if (sub_segments.size() < 4) {
    auto first_seg = sub_segments.begin();
    auto second_seg = sub_segments.begin() + 1;
    // first segment is size 1
    if (sub_segments.begin()->end - sub_segments.begin()->start < min_obs_per_segment) {  // first segment is small
      first_seg = sub_segments.begin();
      second_seg = sub_segments.begin() + 1;
    } else if ((sub_segments.end() - 1)->end - (sub_segments.end() - 1)->start <
               min_obs_per_segment) {  // last segment is small
      first_seg = sub_segments.end() - 2;
      second_seg = sub_segments.end() - 1;
    } else {  // neither edge is "small", can return as-is
      return sub_segments;
    }
    Segment merged_segment = {
        .start = first_seg->start,
        .end = second_seg->end,
        .t = second_seg->t,
        .p = second_seg->p,
        .n = first_seg->n + second_seg->n,
        .mean = second_seg->mean,  // TODO: correct the mean here
    };
    if (sub_segments.size() == 2) {
      res.resize(1);
      res[0] = merged_segment;
      return res;
    } else {  // size == 3
      res.resize(2);
      if (sub_segments.begin()->end - sub_segments.begin()->start < min_obs_per_segment) {  // first segment is small
        res[0] = merged_segment;
        res[1] = sub_segments[2];
      } else if ((sub_segments.end() - 1)->end - (sub_segments.end() - 1)->start <
                 min_obs_per_segment) {  // last segment is small
        res[0] = sub_segments[0];
        res[1] = merged_segment;
      }
      return res;
    }
  } else {
    Logging::Error("Will not prune segment with more than three sub-segments");
    throw std::runtime_error("Will not prune segment with more than three sub-segments");
  }
}
}  // namespace xoos::cnc::segmentation
