#include "segmentation/max-t.h"

#include <tuple>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include "segmentation/cbs-util.h"

namespace xoos::cnc::segmentation {

const std::map<std::string, CbsMaxTMethod> kStringToCbsMethod{
    {"BruteForce", CbsMaxTMethod::kBruteForce},
    {"Fast", CbsMaxTMethod::kFast},
    {"FastMaxN", CbsMaxTMethod::kFastMaxN},
};

/**
 * @brief Gets the sub-segment in obvs that has the maximum t-statistic when compared to the observations outside of the
 * subsegment. This algorithm looks at ALL possible breakpoints (n^2/2) to find the subsegment with the maximum
 * t-statistic.
 * @param obvs observation std::vector
 * @param max_n maximum sub-segment size to test. Set to 0 if no upper limit is needed
 * @param min_n minimum sub-segment size to test.
 * @return Segment contain start and end indices within obvs of the segment w/ that maximizes the t-statistic
 */
Segment GetMaxTBreakpointsBruteForce(const arma::vec& obvs, size_t max_n, size_t min_n) {
  f64 max_t = 0;
  auto m = obvs.n_elem;
  arma::uword max_t_i = 0;
  arma::uword max_t_j = 0;
  arma::vec obvs_psums(arma::cumsum(obvs));
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
      // skip any time either segment is too small
      if (j - i < min_n || i + m - j < min_n) {
        continue;
      }
      f64 t = TStatSquaredMean(obvs_psums, i, j);
      if (std::abs(t) > max_t) {
        max_t = std::abs(t);
        max_t_i = i;
        max_t_j = j;
      }
    }
  }
  Segment ret;
  ret.start = max_t_i;
  ret.end = max_t_j;
  ret.t = max_t;
  return ret;
}

struct PerBlockPartialSumMinMaxRet {
  PerBlockPartialSumMinMaxRet() = default;
  PerBlockPartialSumMinMaxRet(const PerBlockPartialSumMinMaxRet&) = default;
  PerBlockPartialSumMinMaxRet(PerBlockPartialSumMinMaxRet&&) noexcept = default;
  PerBlockPartialSumMinMaxRet& operator=(const PerBlockPartialSumMinMaxRet&) = default;
  PerBlockPartialSumMinMaxRet& operator=(PerBlockPartialSumMinMaxRet&&) noexcept = default;
  arma::uvec block_argmins{};
  arma::uvec block_argmaxs{};
  size_t global_argmin{0};
  size_t global_argmax{0};
};

PerBlockPartialSumMinMaxRet GetPerBlockPartialSumMinMax(const arma::vec& psums, size_t n_blocks, size_t block_size) {
  size_t block_idx = 0;
  arma::uvec block_partial_sum_argmins(n_blocks);
  arma::uvec block_partial_sum_argmaxs(n_blocks);
  block_partial_sum_argmins[0] = 0;
  block_partial_sum_argmaxs[0] = 0;
  f64 block_min_partial_sum = psums[0];
  f64 block_max_partial_sum = psums[0];
  size_t global_min_partial_sum_idx = 0;
  size_t global_max_partial_sum_idx = 0;
  f64 global_min_partial_sum = psums[0];
  f64 global_max_partial_sum = psums[0];
  for (size_t i = 0; i < psums.size(); ++i) {
    if (i != 0 && (i % block_size) == 0) {
      block_idx += 1;
      block_min_partial_sum = psums[i];
      block_max_partial_sum = psums[i];
      block_partial_sum_argmins[block_idx] = i;
      block_partial_sum_argmaxs[block_idx] = i;
    }
    if (psums[i] < block_min_partial_sum) {
      block_partial_sum_argmins[block_idx] = i;
      block_min_partial_sum = psums[i];
    }
    if (psums[i] > block_max_partial_sum) {
      block_partial_sum_argmaxs[block_idx] = i;
      block_max_partial_sum = psums[i];
    }
    if (psums[i] < global_min_partial_sum) {
      global_min_partial_sum_idx = i;
      global_min_partial_sum = psums[i];
    }
    if (psums[i] > global_max_partial_sum) {
      global_max_partial_sum_idx = i;
      global_max_partial_sum = psums[i];
    }
  }
  PerBlockPartialSumMinMaxRet ret;
  ret.block_argmins = std::move(block_partial_sum_argmins);
  ret.block_argmaxs = std::move(block_partial_sum_argmaxs);
  ret.global_argmin = global_min_partial_sum_idx;
  ret.global_argmax = global_max_partial_sum_idx;
  return ret;
}

/**
 * @brief Gets the sub-segment in obvs that has the maximum t-statistic when
 * compared to the observations outside of the subsegment.
 * This algorithm is * greedy, so it reduces the search space of all possible breakpoints at the
 * risk of making the wrong segmentation choice.
 * Algorithm:
 * 1) Divide obvs into sqrt(n) blocks
 * 2) Find the partial sum maxs & mins for each block, and global as well
 * 3) find upper-bound t between all block pairs by comparing each block's partials sum min to the other's partial sum
 * max, then using both options for block boundaries as the lengths of each arc
 * 4) then, consider only the block pairs where the upper-bound t  exceeds max-t
 * 5) update max-t for each block-pair that has a breakpoint pair that exceeds max-t (and update max-t accordingly)
 *
 * @param obvs observation std::vector
 * @param max_n maximum sub-segment size to test. Set to 0 if no upper limit is needed
 * @param min_n minimum sub-segment size to test.
 * @return Segment contain start and end indices within obvs of the segment w/ that maximizes the t-statistic
 */
Segment GetMaxTBreakpointsFast(const arma::vec& obvs, size_t max_n, size_t min_n) {
  arma::vec obvs_psums(arma::cumsum(obvs));
  auto k = static_cast<size_t>(std::ceil(std::sqrt(static_cast<f64>(obvs.size()))));
  PerBlockPartialSumMinMaxRet psum_info = GetPerBlockPartialSumMinMax(obvs_psums, k, k);
  // upper limit for t over the whole array (we call this global)
  auto n = static_cast<f64>(obvs.size());
  size_t max_t_i = std::min(psum_info.global_argmax, psum_info.global_argmin);
  size_t max_t_j = std::max(psum_info.global_argmax, psum_info.global_argmin);
  f64 global_length = std::abs(static_cast<f64>(psum_info.global_argmax) - static_cast<f64>(psum_info.global_argmin));
  f64 global_sum = obvs_psums[psum_info.global_argmax] - obvs_psums[psum_info.global_argmin];
  f64 total_sum = obvs_psums[obvs_psums.size() - 1];
  f64 max_t = TStatSquaredMean(total_sum, n, global_sum, global_length);
  // compare block pairs to find elegible candidates
  using BlockPairKey = std::tuple<size_t, size_t, f64, f64>;
  std::vector<BlockPairKey> candidate_block_pairs;
  // NB: If anything goes wrong with *sensitivity*, then look in this nested for-loop
  for (size_t left_block = 0; left_block < k; ++left_block) {
    size_t left_block_psum_argmin = psum_info.block_argmins[left_block];
    size_t left_block_psum_argmax = psum_info.block_argmaxs[left_block];
    for (size_t right_block = left_block; right_block < k; ++right_block) {
      size_t right_block_psum_argmin = psum_info.block_argmins[right_block];
      size_t right_block_psum_argmax = psum_info.block_argmaxs[right_block];
      // two possible sums, each representing the diff between the max partial sum in one block and the min partial sum
      // in the other
      f64 sum1 = obvs_psums[right_block_psum_argmax] - obvs_psums[left_block_psum_argmin];
      f64 sum2 = obvs_psums[left_block_psum_argmax] - obvs_psums[right_block_psum_argmin];
      // t is probably maximized when the sum differences between the two subsegments are most divergent and the length
      // differences are least divergent
      f64 min_length;
      if (right_block > left_block) {
        // starts at last position of left block, ends at first position of right block
        min_length = static_cast<f64>(((right_block * k) + 1) - ((left_block + 1) * k - 1));
      } else {
        min_length = static_cast<f64>(min_n);
      }
      // starts at first position of left block and ends at last position of right block
      auto max_length =
          static_cast<f64>((right_block + 1) * k > obvs.size() - min_n ? obvs.size() - min_n - (left_block * k)
                                                                       : (right_block + 1) * k - (left_block * k));
      // choose the length and sum that maximizes the t-statistic
      f64 upper_bound_t = std::max({TStatSquaredMean(total_sum, n, sum1, max_length),
                                    TStatSquaredMean(total_sum, n, sum1, min_length),
                                    TStatSquaredMean(total_sum, n, sum2, max_length),
                                    TStatSquaredMean(total_sum, n, sum2, min_length)});
      // if this block pair has potential, start off by recording the t-stat between the max partial sum diff
      // between both blocks
      if (upper_bound_t >= max_t) {
        // calculate t at min and max positions, use larger value
        // subtract 1 from the idxs in order to represent the actual segment being tested
        f64 t1 = TStatSquaredMean(obvs_psums, left_block_psum_argmin + 1, right_block_psum_argmax + 1);
        f64 t2 = TStatSquaredMean(obvs_psums, left_block_psum_argmax + 1, right_block_psum_argmin + 1);
        candidate_block_pairs.emplace_back(left_block, right_block, std::max(t1, t2), upper_bound_t);
      }
    }
  }
  // sort in descending order
  std::sort(candidate_block_pairs.begin(),
            candidate_block_pairs.end(),
            [](const BlockPairKey& left, const BlockPairKey& right) { return std::get<2>(left) > std::get<2>(right); });
  // now go through each block pair and search all breakpoint pairs where first breakpoint is in left block and second
  // breakpoint is in right block we will update max_t accordingly and skip block pairs if needed
  for (const auto& [left_block, right_block, t_at_max_sum, t_upper_lim] : candidate_block_pairs) {
    // if the upper limit for this block pair doesn't even exceed the max_t then we can skip this block pair altogether!
    if (t_upper_lim < max_t) {
      continue;
    }
    // Search in a 1-block radius around the left and right blocks, just to be saf
    size_t left_start_boundary = (left_block - 1) * k > 0 ? (left_block - 1) * k : 0;
    size_t left_end_boundary = (left_block + 2) * k > obvs.size() - min_n ? obvs.size() - min_n : (left_block + 2) * k;
    for (size_t i = left_start_boundary; i < left_end_boundary; ++i) {
      // the right breakpoint is contained in the 1-block flank around right_block
      size_t right_start_boundary;
      if ((right_block - 1) > left_block + 1) {
        right_start_boundary = (right_block - 1) * k + 1;
      } else {
        // if both blocks are the same or adjacent, then we search all breakpoints inside both blocks
        right_start_boundary = i + min_n;
      }
      size_t right_end_boundary = (right_block + 2) * k > obvs.size() ? obvs.size() : (right_block + 2) * k;
      if (right_start_boundary > obvs.size() || right_end_boundary > obvs.size()) {
        throw error::Error("max-t: search boundary exceeds observation vector size");
      }
      auto inner_loop_cond = [right_end_boundary, max_n](size_t i, size_t j) {
        if (max_n > 0) {
          // if max_n is specified, then stop when the segment exceeds this mch
          return j < right_end_boundary && j - i < max_n;
        } else {
          return j < right_end_boundary;
        }
      };
      for (size_t j = right_start_boundary; inner_loop_cond(i, j); ++j) {
        if (i >= obvs.n_elem || j == 0 || j > obvs.n_elem || i >= j) {
          throw error::Error("max-t: breakpoint indices out of bounds");
        }
        // skip any time either segment is too small
        if ((j - i) < min_n || (obvs.size() - j + i) < min_n) {
          continue;
        }
        f64 t = TStatSquaredMean(obvs_psums, i, j);
        if (std::abs(t) > max_t) {
          max_t = std::abs(t);
          max_t_i = i;
          max_t_j = j;
        }
      }
    }
  }
  Segment ret;
  ret.start = max_t_i;
  ret.end = max_t_j;
  ret.t = max_t;
  return ret;
}

/**
 * @brief Like the Fast algorithm but optimized for a small max n
 * @param obvs
 * @param max_n
 * @param
 * @return
 */
Segment GetMaxTBreakpointsFastMaxN(const arma::vec& obvs, size_t max_n, size_t min_n) {
  arma::vec obvs_psums(arma::cumsum(obvs));
  auto n = static_cast<f64>(obvs.size());
  auto n_blocks = static_cast<size_t>(std::ceil(n / static_cast<f64>(max_n)));
  auto psum_info = GetPerBlockPartialSumMinMax(obvs_psums, n_blocks, max_n);
  //  some special logic for comparing segments that cross the last block into the first
  size_t pblock_psum_argmin = psum_info.block_argmins[n_blocks - 1];
  size_t pblock_psum_argmax = psum_info.block_argmaxs[n_blocks - 1];
  f64 pblock_psum_min = obvs_psums[pblock_psum_argmin];
  f64 pblock_psum_max = obvs_psums[pblock_psum_argmax];
  f64 block_psum_min = obvs_psums[0];
  f64 block_psum_max = obvs_psums[0];
  f64 max_t = 0;
  size_t max_t_i = 0;
  size_t max_t_j = 0;
  f64 max_sum_between_blocks =
      std::max(std::abs((obvs_psums[obvs_psums.size() - 1] - pblock_psum_max) + block_psum_min),
               std::abs((obvs_psums[obvs_psums.size() - 1] - pblock_psum_min) + block_psum_max));
  f64 total_sum = obvs_psums[obvs.size() - 1];
  for (size_t length = min_n; length <= max_n; ++length) {
    // calculate the upperbound for the segments crossing the last block and the first block, with this length
    f64 upper_bound_t_between_blocks = TStatSquaredMean(total_sum, obvs.size(), max_sum_between_blocks, length);
    if (upper_bound_t_between_blocks > max_t) {
      for (size_t i = obvs.size() - length + 1; i < obvs.size(); ++i) {
        size_t j = length - (obvs.size() - i);
        if (j + (obvs.size() - i) != length) {
          throw error::Error("max-t: wrap-around segment length mismatch");
        }
        f64 sum_at_first_block = obvs_psums[j - 1];
        // potential bug if first block is equal to last block (is this even possible?)
        f64 sum_at_last_block = obvs_psums[obvs_psums.size() - 1] - (i > 0 ? obvs_psums[i - 1] : 0);
        f64 sum = sum_at_first_block + sum_at_last_block;
        f64 t = TStatSquaredMean(total_sum, obvs.size(), sum, length);
        if (abs(t) > max_t) {
          max_t = abs(t);
          max_t_i = i;
          max_t_j = j;
        }
      }
    }
  }
  for (size_t block = 0; block < n_blocks; ++block) {
    size_t block_psum_argmin = psum_info.block_argmins[block];
    size_t block_psum_argmax = psum_info.block_argmaxs[block];
    block_psum_min = obvs_psums[block_psum_argmin];
    block_psum_max = obvs_psums[block_psum_argmax];
    f64 max_sum_within_block = block_psum_max - block_psum_min;
    f64 max_sum_between_blocks;
    max_sum_between_blocks =
        std::max(std::abs(block_psum_min - pblock_psum_max), std::abs(block_psum_max - pblock_psum_min));
    // for segments inside the block
    for (size_t length = min_n; length <= max_n; ++length) {
      // calculate the within-block upperbound for t for this length
      f64 upper_bound_t_within_block = TStatSquaredMean(total_sum, obvs.size(), max_sum_within_block, length);
      if (upper_bound_t_within_block > max_t) {
        size_t start = block * max_n;
        size_t end = (block + 1) * max_n - length > obvs.size() ? obvs.size() - length : (block + 1) * max_n - length;
        for (size_t i = start; i < end && i + length <= obvs.size(); ++i) {
          f64 t = TStatSquaredMean(obvs_psums, i, i + length);
          if (abs(t) > max_t) {
            max_t = abs(t);
            max_t_i = i;
            max_t_j = i + length;
          }
        }
      }
    }
    // for segments between this block and last
    for (size_t length = min_n; length <= max_n; ++length) {
      // calculate the upperbound for the segments crossing last block and this block
      f64 upper_bound_t_between_blocks = TStatSquaredMean(total_sum, obvs.size(), max_sum_between_blocks, length);
      if (upper_bound_t_between_blocks > max_t) {
        // this should automatically skip the first block-last block pair (we already looked at this above)
        size_t start = (block * max_n > length) ? (block * max_n) - length : 0;
        size_t end = (block * max_n);
        for (size_t i = start; i < end && i + length <= obvs.size(); ++i) {
          f64 t = TStatSquaredMean(obvs_psums, i, i + length);
          if (abs(t) > max_t) {
            max_t = abs(t);
            max_t_i = i;
            max_t_j = i + length;
          }
        }
      }
    }
    pblock_psum_argmin = block_psum_argmin;
    pblock_psum_argmax = block_psum_argmax;
    pblock_psum_min = block_psum_min;
    pblock_psum_max = block_psum_max;
  }

  Segment ret;
  ret.start = max_t_i;
  ret.end = max_t_j;
  ret.t = max_t;
  return ret;
}

Segment GetMaxTBreakpoints(const arma::vec& obvs, size_t max_n, size_t min_n, CbsMaxTMethod method) {
  if (min_n < 1) {
    throw error::Error("max-t: min_n must be at least 1");
  }
  if (min_n > obvs.n_elem) {
    throw error::Error("max-t: min_n ({}) exceeds observation count ({})", min_n, obvs.n_elem);
  }
  if (method == CbsMaxTMethod::kFastMaxN && max_n == 0) {
    throw error::Error("max-t: max_n must be > 0 for kFastMaxN method");
  }
  switch (method) {
    case CbsMaxTMethod::kBruteForce:
      return GetMaxTBreakpointsBruteForce(obvs, max_n, min_n);
    case CbsMaxTMethod::kFast:
      return GetMaxTBreakpointsFast(obvs, max_n, min_n);
    case CbsMaxTMethod::kFastMaxN:
      return GetMaxTBreakpointsFastMaxN(obvs, max_n, min_n);
    default:
      throw std::runtime_error("bad enum");
  }
}
}  // namespace xoos::cnc::segmentation
