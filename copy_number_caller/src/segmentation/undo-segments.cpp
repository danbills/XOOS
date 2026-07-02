#include "segmentation/undo-segments.h"

#include <algorithm>
#include <cassert>

#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "interval-trees.h"

namespace xoos::cnc::segmentation {

constexpr f64 kMedianAbsoluteDeviationScaleFactor = 1.4826;
constexpr size_t kUndoSegmentsMinNumMedians = 2;

/**
 * @brief Calculate the Median Absolute Deviation (MAD) of a vector of observations
 */
static f64 MedianAbsoluteDeviation(const arma::vec& obvs) {
  // median absolute deviation (MAD) is defined as the median of the absolute deviations from the median
  return kMedianAbsoluteDeviationScaleFactor * arma::median(arma::abs(obvs - arma::median(obvs)));
}

/**
 * @brief Calculate the trimmed standard deviation of a vector of observations
 */
f64 GetTrimmedStandardDeviation(const arma::vec& arr) {
  // trimmed.SD <- mad(diff(genomdat))/sqrt(2)
  const arma::vec diff = arr(arma::span(1, arr.n_elem - 1)) - arr(arma::span(0, arr.n_elem - 2));
  f64 trimmed_sd = MedianAbsoluteDeviation(diff);
  trimmed_sd = trimmed_sd / std::pow(2.0, 0.5);
  return trimmed_sd;
}

/**
 * @brief holds segment with its associated observations, for each of iteration
 */
struct SegmentWithObservations {
  GenomicSegment segment;
  Observations obvs;
};

Observations MergeObservations(const Observations& obvs1, const Observations& obvs2) {
  Observations merged_obvs;
  merged_obvs.contigs.reserve(obvs1.contigs.size() + obvs2.contigs.size());
  merged_obvs.contigs.insert(merged_obvs.contigs.begin(), obvs1.contigs.begin(), obvs1.contigs.end());
  merged_obvs.contigs.insert(merged_obvs.contigs.end(), obvs2.contigs.begin(), obvs2.contigs.end());
  merged_obvs.starts = arma::join_cols(obvs1.starts, obvs2.starts);
  merged_obvs.ends = arma::join_cols(obvs1.ends, obvs2.ends);
  merged_obvs.obvs = arma::join_cols(obvs1.obvs, obvs2.obvs);
  return merged_obvs;
}

/**
 * @brief Undo genomic segments by merging adjacent segments if their median observation difference is less than a
 * factor (max_trimmed_sd_factor) of the trimmed standard deviation (trimmed_sd)
 * @param segments input segments with associated observations
 * @param trimmed_sd trimmed standard deviation of the observations (use GetTrimmedStandardDeviation to calculate)
 * @param max_trimmed_sd_factor maximum factor of trimmed standard deviation to consider for merging segments
 * @return merged segments
 */
static std::vector<SegmentWithObservations> UndoGenomicSegments(const std::vector<SegmentWithObservations>& segments,
                                                                f64 trimmed_sd,
                                                                f64 max_trimmed_sd_factor) {
  if (segments.empty()) {
    return {};
  }
  f64 max_threshold = trimmed_sd * max_trimmed_sd_factor;
  // merging will be done in a series of iterations. In each iteration, we will transfer segments from
  // `prev_merged_segments` to `merged_segments`, merging segments in `prev_merged_segments` if their median observation
  // difference is less than the threshold, and then compare the size of `merged_segments` with that of
  // `prev_merged_segments` to determine whether to continue iterating (i.e. whether any merging has been done in this
  // iteration)
  std::vector<SegmentWithObservations> prev_merged_segments(segments);
  std::vector<SegmentWithObservations> merged_segments;
  merged_segments.reserve(prev_merged_segments.size());
  std::vector<f64> medians;
  medians.reserve(segments.size());
  std::vector<u8> breakpoints_to_merge(segments.size());
  while (merged_segments.size() != prev_merged_segments.size()) {
    merged_segments.clear();
    medians.clear();
    // calculate the medians of all the segments
    for (const auto& seg : prev_merged_segments) {
      medians.emplace_back(arma::median(seg.obvs.obvs));
    }
    // this vector will tell us which segments to merge together.  0 means don't merge, 1 means merge
    breakpoints_to_merge.clear();
    breakpoints_to_merge.resize(medians.size());
    // breakpoints[i] == 1 will imply that we need to merge semgnet[i] with segment[i-1]. breakpoints[i] == 0 means do
    // not merge segment[i] with segment[i-1].
    std::ranges::fill(breakpoints_to_merge.begin(), breakpoints_to_merge.end(), 0);
    // for each element in medians, compare to the previous and next element
    for (size_t i = 1; i < medians.size() - 1; ++i) {
      const f64 prev_diff = std::abs(medians[i] - medians[i - 1]);
      const f64 next_diff = std::abs(medians[i] - medians[i + 1]);
      // make sure the distance to this neighbor is less than the difference to the next neighbor. Cannot merge both
      // neighbors in same pass
      if (prev_diff <= next_diff && prev_diff <= max_threshold) {
        // if it is 1, leave it as 1. If it is 0, change to 1
        breakpoints_to_merge[i] += 1 - breakpoints_to_merge[i];
      } else if (next_diff < prev_diff && next_diff <= max_threshold) {
        breakpoints_to_merge[i + 1] += 1 - breakpoints_to_merge[i + 1];
      } else {
        // if neither neighbor is less than the threshold, then we don't merge this segment with either neighbor, so do
        // nothing and leave the value as is
      }
    }
    if (medians.size() >= kUndoSegmentsMinNumMedians) {
      // handle the last breakpoint
      const f64 prev_diff = std::abs(medians[medians.size() - 1] - medians[medians.size() - 2]);
      if (prev_diff <= max_threshold) {
        // if it is 1, leave it as 1. If it is 0, change to 1
        breakpoints_to_merge[medians.size() - 1] += 1 - breakpoints_to_merge[medians.size() - 1];
      }
    }
    // begin merge process. We iteratively merge segments with their left
    // first segment automatically added to merged_segments because it can't be merged to the left.
    merged_segments.emplace_back(prev_merged_segments[0]);
    // walk through the breakpoints and medians, and simultaneously merge segments from prev_merged_segments and add to
    // the back of merged_segments if breakpoint is 1, otherwise just add the segment to merged_segments without merging
    for (size_t i = 1; i < medians.size(); ++i) {
      if (breakpoints_to_merge[i] > 0) {
        merged_segments.back().segment.end = prev_merged_segments[i].segment.end;
        // also transfer all the observations of the to-be-merged segment
        merged_segments.back().obvs = MergeObservations(merged_segments.back().obvs, prev_merged_segments[i].obvs);
      } else {
        merged_segments.emplace_back(prev_merged_segments[i]);
      }
    }
    // in the next round, we will merge the segments even further (if possible), so swap merged_segments with
    // prev_merged_segments. merged_segments will be cleared at the beginning of the next loop
    std::swap(prev_merged_segments, merged_segments);
  }
  return prev_merged_segments;
}

/**
 * @brief wrapper around UndoSegments to work with GenomicSegment and Observations
 */
std::vector<GenomicSegment> UndoSegments(const std::vector<GenomicSegment>& segments,
                                         const Observations& obvs,
                                         const f64 trimmed_sd,
                                         const f64 max_trimmed_sd_factor) {
  // convert segments and obvs into SegmentWithObservations
  std::vector<SegmentWithObservations> segs_with_obvs;
  const IntervalTrees obvs_trees(obvs);
  for (const auto& seg : segments) {
    const Observations seg_obvs = obvs.FilterByIdxs(obvs_trees.LookUp(seg.contig, seg.start, seg.end));
    segs_with_obvs.emplace_back(SegmentWithObservations{seg, seg_obvs});
  }
  // call UndoGenomicSegments
  const std::vector<SegmentWithObservations> undone_segs =
      UndoGenomicSegments(segs_with_obvs, trimmed_sd, max_trimmed_sd_factor);
  // extract GenomicSegments from undone_segs
  std::vector<GenomicSegment> ret;
  ret.reserve(undone_segs.size());
  for (const auto& seg_with_obvs : undone_segs) {
    ret.emplace_back(seg_with_obvs.segment);
  }
  return ret;
}

/**
 * @brief Iteratively undo segments until the number of segments is less than or equal to max_num_segments. Increase
 * max_trimmed_sd_factor by increment_max_trimmed_sd_factor at each iteration. Note that except for the contig, start,
 * and end fields of each GenomicSegment memebr, all other fields in these members are invalidated
 * @param seed segments seed segments
 * @param segments input segments
 * @param obvs observations vector
 * @param trimmed_sd trimmed standard deviation of the observations
 * @param start_max_trimmed_sd_factor starting max trimmed sd factor
 * @param increment_max_trimmed_sd_factor increment max trimmed sd factor at each iteration
 * @param max_num_segments maximum number of segments to achieve
 * @return merged segments
 */
std::vector<GenomicSegment> IterativeUndoGenomicSegments(const std::vector<GenomicSegment>& seed_segments,
                                                         const std::vector<GenomicSegment>& segments,
                                                         const Observations& obvs,
                                                         const f64 start_max_trimmed_sd_factor,
                                                         const f64 increment_max_trimmed_sd_factor,
                                                         const size_t max_num_segments) {
  if (segments.size() <= max_num_segments) {
    Logging::Info("No need for iterative undoing: {} <= {} segments", segments.size(), max_num_segments);
    return segments;
  } else {
    Logging::Info("Iteratively undoing segments. Targeting {} -> {} segments", segments.size(), max_num_segments);
  }
  // associate segments with seed segments, and observations with segments
  const IntervalTrees segment_trees(segments);
  const IntervalTrees obvs_trees(obvs);
  std::vector<std::vector<SegmentWithObservations>> segments_per_seed(seed_segments.size());
  for (size_t i = 0; i < seed_segments.size(); ++i) {
    auto segment_idxs = segment_trees.LookUp(seed_segments[i].contig, seed_segments[i].start, seed_segments[i].end);
    for (auto segment_idx : segment_idxs) {
      const GenomicSegment& seg = segments[segment_idx];
      segments_per_seed[i].emplace_back(seg, obvs.FilterByIdxs(obvs_trees.LookUp(seg.contig, seg.start, seg.end)));
    }
  }
  // prepare for undoing segments
  std::vector<SegmentWithObservations> all_undone_segments;
  for (const auto& vec : segments_per_seed) {
    all_undone_segments.insert(all_undone_segments.end(), vec.begin(), vec.end());
  }
  f64 max_trimmed_sd_factor = start_max_trimmed_sd_factor;
  // this iter determines how much to increment max_trimmed_sd_factor at each iteration
  size_t iter = 0;
  size_t prev_num_segments = all_undone_segments.size() + 1;
  // intiate loop; call UndoGenomicSegments on the original segments until new_size <= max_num_segments
  const f64 trimmed_sd = GetTrimmedStandardDeviation(obvs.obvs);
  // get the trimmed standard deviation for each the segments in seg_vec
  for (iter = 0; all_undone_segments.size() > max_num_segments && all_undone_segments.size() < prev_num_segments;
       ++iter) {
    prev_num_segments = all_undone_segments.size();
    all_undone_segments.clear();
    max_trimmed_sd_factor += static_cast<f64>(iter) * increment_max_trimmed_sd_factor;
    for (const auto& seg_vec : segments_per_seed) {
      if (!seg_vec.empty()) {
        std::vector<SegmentWithObservations> undone_segments =
            UndoGenomicSegments(seg_vec, trimmed_sd, max_trimmed_sd_factor);
        all_undone_segments.insert(all_undone_segments.end(), undone_segments.begin(), undone_segments.end());
      }
    }
    Logging::Info("Iteration {}: max_trimmed_sd_factor = {}, number of segments = {}",
                  iter,
                  max_trimmed_sd_factor,
                  all_undone_segments.size());
  }
  Logging::Info("Iterative undoing completed in {} iterations: {}->{} segments",
                iter,
                segments.size(),
                all_undone_segments.size());
  std::vector<GenomicSegment> ret;
  ret.reserve(max_num_segments);
  for (const auto& seg_with_obvs : all_undone_segments) {
    ret.emplace_back(seg_with_obvs.segment);
  }
  return ret;
}
}  // namespace xoos::cnc::segmentation
