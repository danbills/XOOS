/**
 * reference:
 * Olshen et al, 2004. "Circular binary segmentation for the analysis of array‐based DNA copy number data"
 * https://academic.oup.com/biostatistics/article/5/4/557/275197
 */

#include "segmentation/cbs-for-somatic.h"

#include <cmath>
#include <limits>
#include <ranges>  // NOLINT cpplint doesn't yet support some C++20 headers
#include <stack>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/stats/external/asa243/asa243.hpp>
#include <xoos/util/math.h>

#include "segmentation/boundaries.h"

namespace xoos::cnc::segmentation::somatic {

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::chrono::nanoseconds;
using std::chrono::seconds;

const size_t kCbsSiegmundPIters = 100;
const size_t kCbsMinObvsN = 10;
const size_t kCbsHybridK = 25;
const size_t kCbsMaxStackDepth = 1024;

/**
 * @brief Imagine segment as circle joined at both ends. This function calculates t-statistic for testing the hypothesis
 * that arc from i+1 to j and its complement have different means.
 * @param n_v1 size of minor arc
 * @param mean_v1 mean of minor arc
 * @param var_v1 variance of minor arc
 * @param n_v2 size of major arc
 * @param mean_v2 mean of major arc
 * @param var_v2 variance of major arc
 * @return t-statistic describing difference in means between minor and major arc
 */
f64 TStat(f64 n_v1, f64 mean_v1, f64 var_v1, f64 n_v2, f64 mean_v2, f64 var_v2) {
  // pooled variance of v1 and v2
  f64 s_ij = std::pow(((n_v1 - 1) * var_v1 + (n_v2 - 1) * var_v2) / (n_v1 + n_v2 - 2), 0.5);
  // t stat is a function of the means of both std::vectors and their pooled variance
  f64 t_ij = (mean_v1 - mean_v2) / (s_ij * std::pow((1.0 / n_v1 + 1.0 / n_v2), 0.5));
  return t_ij;
}

/**
 * @brief calculates t-statistic between two std::vectors that can be of unequal size and variance
 * @param v1 std::vector 1
 * @param v2 std::vector 2
 * @return t-statistic between v1 and v2
 */
f64 TStatFromVecs(const arma::vec& v1, const arma::vec& v2) {
  auto n_v1 = static_cast<f64>(v1.n_elem);
  auto n_v2 = static_cast<f64>(v2.n_elem);
  f64 mean_v1 = arma::mean(v1);
  f64 mean_v2 = arma::mean(v2);
  f64 var_v1 = arma::var(v1);
  f64 var_v2 = arma::var(v2);
  return TStat(n_v1, mean_v1, var_v1, n_v2, mean_v2, var_v2);
}

/* warning - this is the biased estimator for variance. Do we need an unbiased estimator??? */
static inline f64 Variance(f64 m, f64 mean, f64 sq_sum) {
  return (sq_sum / m) - std::pow(mean, 2);
}

/**
 * @brief Gets the sub-segment in obvs that has the maximum t-statistic when compared to the observations outside of the
 * subsegment
 * @param obvs observation std::vector
 * @param max_n maximum sub-segment size to test. Set to 0 if no upper limit is needed
 * @param min_n minimum sub-segment size to test.
 * @return Segment contain start and end indices within obvs of the segment w/ that maximizes the t-statistic
 */
Segment GetMaxTSegCumSum(const arma::vec& obvs, size_t max_n, size_t min_n) {
  f64 max_t = 0;
  f64 max_t_mean = 0;
  f64 max_t_m = 0;
  auto m = obvs.n_elem;
  arma::uword max_t_i = 0;
  arma::uword max_t_j = 0;
  arma::vec obvs_psums(arma::cumsum(obvs));
  arma::vec obvs_psqsums(arma::cumsum(arma::pow(obvs, 2)));
  for (arma::uword i = 0; i < obvs.n_elem - 1; ++i) {
    auto inner_loop_cond = [&obvs, max_n](arma::uword i, arma::uword j) {
      if (max_n > 0) {
        // if max_n is specified, then stop when the segment exceeds this mch
        return j < obvs.n_elem && j - i < max_n;
      } else {
        return j < obvs.n_elem;
      }
    };
    for (arma::uword j = i + min_n; inner_loop_cond(i, j); ++j) {
      auto m_v1 = static_cast<f64>(j - i);
      auto m_v2 = static_cast<f64>(i + m - j);
      // skip any time either segment is too small
      if (j - i < min_n || i + m - j < min_n) {
        continue;
      }
      f64 mean_v1 = (obvs_psums[j - 1] - (i > 0 ? obvs_psums[i - 1] : 0)) / m_v1;
      f64 var_v1 = Variance(m_v1, mean_v1, obvs_psqsums[j - 1] - (i > 0 ? obvs_psqsums[i - 1] : 0));
      f64 mean_v2;
      f64 var_v2;
      if (i > 0 && j < obvs.n_elem) {
        mean_v2 = (obvs_psums[i - 1] + obvs_psums[m - 1] - obvs_psums[j - 1]) / m_v2;
        var_v2 = Variance(m_v2, mean_v2, obvs_psqsums[i - 1] + obvs_psqsums[m - 1] - obvs_psqsums[j - 1]);
      } else if (i == 0 && j < m) {
        mean_v2 = (obvs_psums[m - 1] - obvs_psums[j - 1]) / m_v2;
        var_v2 = Variance(m_v2, mean_v2, obvs_psqsums[m - 1] - obvs_psqsums[j - 1]);
      } else if (i > 0 && j == m) {
        mean_v2 = (obvs_psums[i - 1]) / m_v2;
        var_v2 = Variance(m_v2, mean_v2, obvs_psqsums[i - 1]);
      } else {
        throw std::runtime_error("TStat: i==0 and j==m; cannot make segment");
      }
      f64 t = TStat(m_v1, mean_v1, var_v1, m_v2, mean_v2, var_v2);
      if (std::abs(t) > max_t) {
        max_t = std::abs(t);
        max_t_i = i;
        max_t_j = j;
        max_t_mean = mean_v1;
        max_t_m = m_v1;
      }
    }
  }
  Segment ret;

  ret.start = max_t_i;
  ret.end = max_t_j;
  ret.t = max_t;
  ret.mean = max_t_mean;
  ret.n = max_t_m;
  return ret;
}

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
 * @param min_obs_size minimum number of observations allowed per segment
 * @return true if null is rejected else false
 */
bool CircularPermuteNullAndTestSignificance(const arma::vec& obvs,
                                            size_t n_permutations,
                                            f64 max_p,
                                            f64 t,
                                            size_t max_obs_size,
                                            const arma::uvec& boundaries,
                                            size_t min_obs_size) {
  auto max_n_gt_t = static_cast<size_t>(max_p * static_cast<f64>(n_permutations));
  Logging::Debug("Estimating null over maximum {} iters using early stopping measure. hard stop at {} excds",
                 n_permutations,
                 max_n_gt_t + 1);
  size_t n_gt_t = 0;
  // k denotes start-index in boundaries std::vector
  size_t k = max_n_gt_t * (max_n_gt_t + 1) / 2;
  Logging::Debug("boundary (std::vector size {}) start index is: {}", boundaries.n_elem, k);
  for (size_t i = 0; i < n_permutations; ++i) {
    arma::arma_rng::set_seed(i);
    arma::vec obvs_perm = arma::shuffle(obvs);  // is there a faster way to iterate through permutations?
    Segment perm_max_t_seg = GetMaxTSegCumSum(obvs_perm, max_obs_size, min_obs_size);
    if (t < perm_max_t_seg.t) {
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
bool BinaryPermuteNullAndTestSignificance(
    const arma::vec& obvs, size_t n_permutations, f64 max_p, const Segment& seg1, const Segment& seg2) {
  arma::vec to_cmp =
      arma::join_cols(obvs(arma::span(seg1.start, seg1.end - 1)), obvs(arma::span(seg2.start, seg2.end - 1)));
  arma::vec v1 = obvs(arma::span(seg1.start, seg1.end - 1));
  arma::vec v2 = obvs(arma::span(seg2.start, seg2.end - 1));
  size_t n_v1 = v1.n_elem;
  f64 nexcd = 0.0;
  f64 max_excd = static_cast<f64>(n_permutations) * max_p + 1;
  f64 t = std::abs(TStatFromVecs(v1, v2));
  for (size_t i = 0; i < n_permutations; ++i) {
    arma::vec perm = arma::shuffle(to_cmp);
    v1 = perm(arma::span(0, n_v1 - 1));
    v2 = perm(arma::span(n_v1, perm.n_elem - 1));
    f64 perm_t = std::abs(TStatFromVecs(v1, v2));
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
    throw error::Error("cbs-for-somatic: SiegmundP trapz returned unexpected number of elements");
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
 * @param min_obs_size minimum number of observations allowed per segment
 * @return true if the segmentation was significant according to the p-value calculation, false if not
 */
bool HybridTestForSignificance(
    const arma::vec& obvs, size_t n_permutations, f64 max_p, f64 t, const arma::uvec& boundaries, size_t min_obs_size) {
  if (obvs.n_elem < kCbsHybridK * 4) {
    Logging::Debug("{} elem segment: using default P method", obvs.n_elem);
    return CircularPermuteNullAndTestSignificance(obvs, n_permutations, max_p, t, 0, boundaries, min_obs_size);
  } else {
    Logging::Debug("{} elem segment: using approx P method", obvs.n_elem);
    // get null distribution for all splits minor_arc>=k && major_arc>=k
    f64 p = SiegmundP(t, static_cast<f64>(obvs.n_elem), static_cast<f64>(kCbsHybridK));
    if (p < max_p) {
      // get null distribution for all splits minor_arc<k to verify that p actually is significant
      Logging::Debug("{} elem segment: approx p is significant ({}), will also run P method", obvs.n_elem, p);
      return CircularPermuteNullAndTestSignificance(
          obvs, n_permutations, max_p - p, t, kCbsHybridK, boundaries, min_obs_size);
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
 * @param min_obs_size minimum number of observations allowed per segment
 * @return a maximum of three segments w.r.t to obvs. These segments should span the entirety of obvs
 */
std::vector<Segment> GetSegments(
    const arma::vec& obvs, f64 max_p, size_t n_permutations, const arma::uvec& boundaries, size_t min_obs_size) {
  // get the highest t-stat and corresponding segment for the real observations
  Logging::Debug("finding max-T sub-segment for {} obvs segment", obvs.n_elem);
  Segment max_t_seg = GetMaxTSegCumSum(obvs, 0, min_obs_size);
  Logging::Debug("found max-T sub-segment for {} obvs segment", obvs.n_elem);
  if (math::IsCloseToZero(max_t_seg.t)) {
    Logging::Debug("segment with {} obvs is completely uniform", obvs.n_elem);
    return {};
  }
  if (max_t_seg.end - max_t_seg.start < kCbsMinObvsN) {
    Logging::Debug("segment with {} obvs has very small max-t sub-segment", obvs.n_elem);
    return {};
  }
  bool is_significant = false;
  is_significant = HybridTestForSignificance(obvs, n_permutations, max_p, max_t_seg.t, boundaries, min_obs_size);
  Logging::Debug("found significance for max-T sub-segment for {} obvs segment", obvs.n_elem);
  if (!is_significant) {
    return std::vector<Segment>{};
  }
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
 * @param min_obs_size minimum number of observations allowed per segment
 * @return all significant segments found in obvs
 */
std::vector<Segment> CircularBinarySegmentation(
    const arma::vec& obvs, f64 max_p, size_t n_permutations, size_t min_obs_size, bool no_single_obvs_subsegments) {
  auto time1 = high_resolution_clock::now();
  if (obvs.empty()) {
    return {};
  }
  std::vector<Segment> ret;
  std::stack<Segment> seg_stack;
  Segment first_seg = {.start = 0, .end = obvs.n_elem};
  seg_stack.push(first_seg);
  while (!seg_stack.empty()) {
    Logging::Debug("stack has {} elements", seg_stack.size());
    const Segment& segment = seg_stack.top();
    seg_stack.pop();
    // don't process a segment if it's too small
    if (segment.end - segment.start < kCbsMinObvsN) {
      Logging::Debug("adding small segment with {} obvs", (segment.end - segment.start));
      ret.push_back(segment);
    } else {
      arma::vec sub_obvs = obvs(arma::span(segment.start, segment.end - 1));
      // "recursive" step - partition segment into sub-segments
      arma::uvec boundaries;
      if (n_permutations == 10000 && math::IsEqualWithTolerance(max_p, 0.01)) {
        boundaries = kBoundaryP10000A01Eta05;
      }
      Logging::Debug("Segmenting on {} obs", sub_obvs.n_elem);
      std::vector<Segment> sub_segments = GetSegments(sub_obvs, max_p, n_permutations, boundaries, min_obs_size);
      Logging::Debug("Segmented on {} obs", sub_obvs.n_elem);
      // merge back single-obs sub-segments
      if (no_single_obvs_subsegments) {
        sub_segments = MergeSingleObvsEndSegment(sub_segments, min_obs_size);
      }
      if (sub_segments.size() <= 1) {
        // if this segment cannot be broken down, then return it as one segment
        Logging::Debug("adding un-breakable segment with {} obvs", (segment.end - segment.start));
        ret.push_back(segment);
      } else {
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
          if (!BinaryPermuteNullAndTestSignificance(sub_obvs, n_permutations, max_p, seg_to_merge1, seg_to_merge2)) {
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
            if (new_segment.end > obvs.n_elem) {
              throw error::Error(
                  "cbs-for-somatic: sub-segment end exceeds observation count after coordinate adjustment");
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
  }
  std::sort(ret.begin(), ret.end(), [](const Segment& l, const Segment& r) { return l.start < r.start; });
  auto time2 = high_resolution_clock::now();
  auto total_time = duration_cast<milliseconds>(time2 - time1).count();
  Logging::Debug("found {} segments over {} elements in {}ms", ret.size(), obvs.n_elem, total_time);
  return ret;
}

/**
 * @brief Given a set of 2-3 subsegments, if the first or the last sub-segment has fewer than min_obs_size observations,
 * then merge it with its adjacent sub-segment
 * @param sub_segments a std::vector of 2-3 sub-segments
 * @param min_obs_size  minimum number of observations allowed per segment
 * @return
 */
std::vector<Segment> MergeSingleObvsEndSegment(const std::vector<Segment>& sub_segments, size_t min_obs_size) {
  if (sub_segments.size() <= 1) {
    return sub_segments;
  }
  std::vector<Segment> res;
  if (sub_segments.size() < 4) {
    auto first_seg = sub_segments.begin();
    auto second_seg = sub_segments.begin() + 1;
    // first segment is size 1
    if (sub_segments.begin()->end - sub_segments.begin()->start < min_obs_size) {  // first segment is small
      first_seg = sub_segments.begin();
      second_seg = sub_segments.begin() + 1;
    } else if ((sub_segments.end() - 1)->end - (sub_segments.end() - 1)->start <
               min_obs_size) {  // last segment is small
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
      if (sub_segments.begin()->end - sub_segments.begin()->start < min_obs_size) {  // first segment is small
        res[0] = merged_segment;
        res[1] = sub_segments[2];
      } else if ((sub_segments.end() - 1)->end - (sub_segments.end() - 1)->start <
                 min_obs_size) {  // last segment is small
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
}  // namespace xoos::cnc::segmentation::somatic
