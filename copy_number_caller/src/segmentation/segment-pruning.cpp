#include "segmentation/segment-pruning.h"

#include <xoos/log/logging.h>

#include "segmentation/fastcluster-helper.h"
#include "segmentation/interval-trees.h"

namespace xoos::cnc::segmentation {
/**
 * @brief breaks a set of logr genomicsegmnets according to the partitioning defined by dh_segments. for a
 set of dh segments overlapping a logr segment, the dh segments are extended such that the endpoints of each segment
 correspond to the LogR observation that surrounds the midpoint between the two BAF segments.
 * @param logr_segments
 * @param logrs
 * @param dh_segments
 * @param dhs
 * @return a list of genomicsegment the same size as dh_segments, with the boundaries extended to fully cover the
 genomic space of logr_segments

 */
std::vector<GenomicSegment> BreakLogrSegmentsByDHSegments(const std::vector<GenomicSegment>& logr_segments,
                                                          const Observations& logrs,
                                                          const std::vector<GenomicSegment>& dh_segments,
                                                          const Observations& dhs) {
  std::vector<GenomicSegment> new_segments = dh_segments;
  IntervalTrees dh_trees(dh_segments);
  for (const auto& logr_seg : logr_segments) {
    // find all the dh segments that overlap this logr segment
    std::vector<size_t> dh_seg_idxs = dh_trees.LookUp(logr_seg.contig, logr_seg.start, logr_seg.end);
    // special handling when there are no or just one BAF segment overlapping logr_seg
    if (dh_seg_idxs.size() == 1) {
      new_segments[dh_seg_idxs[0]].start = logr_seg.start;
      new_segments[dh_seg_idxs[0]].end = logr_seg.end;
      continue;
    } else if (dh_seg_idxs.empty()) {
      continue;
    }
    new_segments[dh_seg_idxs[0]].start = logr_seg.start;
    size_t left_idx = 0;
    size_t right_idx = 1;
    size_t left_seg_idx = dh_seg_idxs[left_idx];
    size_t right_seg_idx = dh_seg_idxs[right_idx];
    while (right_idx < dh_seg_idxs.size()) {
      // genomic position between the end of the left DH segment and the start of the right DH segment
      size_t midpoint = (new_segments[left_seg_idx].end + new_segments[right_seg_idx].start) / 2;
      // find the logr value in logr_seg that precedes midpoint
      for (size_t logr_idx = logr_seg.arr_start.value() + 1; logr_idx < logr_seg.arr_end.value(); ++logr_idx) {
        if (logrs.ends[logr_idx] >= midpoint) {
          // If it's the case that the left segment falls entirely within the
          // observation, keep merging it with the right segment until the right
          // end is greater than the end of the observation
          if (new_segments[left_seg_idx].start >= logrs.starts[logr_idx]) {
            while (midpoint < logrs.ends[logr_idx] && right_idx < dh_seg_idxs.size()) {
              midpoint = (new_segments[left_seg_idx].end + new_segments[right_seg_idx].start) / 2;
              new_segments[left_seg_idx].end = new_segments[right_seg_idx].end;
              new_segments[left_seg_idx].arr_end = new_segments[right_seg_idx].arr_end;
              new_segments[right_seg_idx].start = 0;
              new_segments[right_seg_idx].end = 0;
              new_segments[right_seg_idx].arr_start = 0;
              new_segments[right_seg_idx].arr_end = 0;
              right_idx++;
              right_seg_idx = dh_seg_idxs[right_idx];
            }
            if (right_idx == dh_seg_idxs.size()) {
              break;
            }
            // since now the midpoint is (almost) "past" the current logr segment, we need to increment logr_idx
            logr_idx += 1;
            // if it's the last logr idx, then we extend the left segment to the end, invalidate the right segment, and
            // stop
            if (logr_idx == logr_seg.arr_end.value()) {
              new_segments[left_seg_idx].end = logrs.ends[logr_seg.arr_end.value() - 1];
              new_segments[left_seg_idx].arr_end = new_segments[right_seg_idx].arr_end;
              new_segments[right_seg_idx].start = 0;
              new_segments[right_seg_idx].end = 0;
              break;
            }
          }
          // end of left segment should be end of PREVIOUS logr observation
          new_segments[left_seg_idx].end = logrs.ends[logr_idx - 1];
          // start of right segment should be start of this logr observation
          new_segments[right_seg_idx].start = logrs.starts[logr_idx];
          break;
        }
      }
      left_idx = right_idx;
      left_seg_idx = dh_seg_idxs[left_idx];
      right_idx = left_idx + 1;
      right_seg_idx = dh_seg_idxs[right_idx];
    }
  }
  std::erase_if(new_segments, [](const GenomicSegment& gs) { return gs.start == 0 && gs.end == 0; });
  // re-scan new segments and make sure that logr-segments are _fully covered_. If not, extend the segments to the logr
  // segments boundaries
  IntervalTrees new_segment_trees(new_segments);
  for (const auto& logr_seg : logr_segments) {
    std::vector<size_t> new_seg_idxs = new_segment_trees.LookUp(logr_seg.contig, logr_seg.start, logr_seg.end);
    // inspect first and last new segments. Adjust their breakpoints if necessary
    if (!new_seg_idxs.empty()) {
      if (new_segments[new_seg_idxs.front()].start > logr_seg.start) {
        new_segments[new_seg_idxs.front()].start = logr_seg.start;
      }
      if (new_segments[new_seg_idxs.back()].end < logr_seg.end) {
        new_segments[new_seg_idxs.back()].end = logr_seg.end;
      }
    }
  }
  return new_segments;
}

/**
 * @brief Pruning segments by clustering segments with similar mean logr and weighted SNP AF if provided
 * @param logr_segments
 * @param joint_clustering bool, whether or not to perform join clustering
 * @param h height parameter in Hierarchical clustering
 * @return a list of pruned genomicsegment
 */
std::vector<GenomicSegment> PruneSegmentsByClustering(const std::vector<GenomicSegment>& logr_segments,
                                                      const bool joint_clustering,
                                                      std::optional<f64> h) {
  auto time1 = std::chrono::high_resolution_clock::now();
  constexpr s32 kMinSnps = 5;
  std::vector<GenomicSegment> new_segments;

  size_t n_segments = logr_segments.size();
  if (n_segments <= 1) {
    Logging::Info("Less than two segments found, skip Pruning...", n_segments);
    new_segments = logr_segments;
    return new_segments;
  } else {
    Logging::Info("Pruning {} segments with hierarchical clustering", n_segments);
  }
  std::vector<bool> ignore_seg(n_segments, false);
  std::vector<Point> points;

  for (size_t k = 0; k < n_segments; k++) {
    ignore_seg[k] = !logr_segments[k].num_obs.has_value() || logr_segments[k].num_obs == 0 ||
                    !logr_segments[k].mean_logr.has_value() || std::isnan(logr_segments[k].mean_logr.value());
    if (joint_clustering) {
      ignore_seg[k] = ignore_seg[k] && (!logr_segments[k].num_snps.has_value() || logr_segments[k].num_snps < kMinSnps);
    }
  }

  // Some segments have no logr, ignore these segments for clustering, but maintain an index to retrieve them
  std::vector<s32> original_idx(n_segments, -1);
  s32 k = 0;
  f64 mean_af = 0.0;
  for (size_t i = 0; i < n_segments; i++) {
    if (logr_segments[i].mean_logr.has_value() && !std::isnan(logr_segments[i].mean_logr.value())) {
      if (joint_clustering && logr_segments[i].mean_dh.has_value() && !std::isnan(logr_segments[i].mean_dh.value())) {
        mean_af = logr_segments[i].mean_dh.value() / 2.0;  // mirrored AF = 1/2 DH
      }
      points.emplace_back(logr_segments[i].mean_logr.value(), mean_af);
      original_idx[i] = k++;
    }
  }

  // TODO: infer h based on data
  f64 h_new = 0.25;
  if (h.has_value()) {
    h_new = h.value();
  }
  Logging::Debug("Setting clustering paramter H = {}", h_new);

  std::vector<s32> clusters = Hclust(points, hclust_fast_methods::HCLUST_METHOD_COMPLETE, h_new);

  Logging::Debug("Merging segments based on clustering");

  // Setting segment mean logr as weighted mean of segments in the same cluster, weighted by num of snps
  std::vector<f64> logr_mean(clusters.size(), 0.0);
  std::vector<s32> total_num_snp(clusters.size(), 0);
  std::vector<f64> logr_mean_sum(clusters.size(), 0.0);
  for (size_t i = 0; i < n_segments; i++) {
    if (ignore_seg[i] || original_idx[i] < 0) {
      continue;
    }
    s32 cluster_id = clusters[original_idx[i]];
    total_num_snp[cluster_id] += static_cast<s32>(logr_segments[i].num_obs.value());
    logr_mean_sum[cluster_id] +=
        logr_segments[i].mean_logr.value() * static_cast<s32>(logr_segments[i].num_obs.value());
  }

  for (size_t i = 0; i < total_num_snp.size(); i++) {
    if (total_num_snp[i] > 0) {
      logr_mean[i] = logr_mean_sum[i] / total_num_snp[i];
    }
  }

  // Mering consecutive segments with the same cluster id
  GenomicSegment gs = logr_segments[0];
  for (size_t i = 1; i < n_segments; i++) {
    if (ignore_seg[i - 1] || ignore_seg[i] || logr_segments[i].contig != logr_segments[i - 1].contig ||
        original_idx[i] < 0 || original_idx[i - 1] < 0 || clusters[original_idx[i]] != clusters[original_idx[i - 1]]) {
      new_segments.emplace_back(gs);
      gs = logr_segments[i];
    } else {
      gs.end = logr_segments[i].end;
      gs.arr_end = logr_segments[i].arr_end;
      gs.num_obs = gs.num_obs.value() + logr_segments[i].num_obs.value();
      gs.mean_logr = logr_mean[clusters[original_idx[i]]];
      if (gs.num_snps.has_value() && logr_segments[i].num_snps.has_value()) {
        gs.num_snps = gs.num_snps.value() + logr_segments[i].num_snps.value();
      }
      // mean_dh is invalidated by merging; clear it so callers re-populate from raw observations
      gs.mean_dh = std::nullopt;
    }
  }
  new_segments.emplace_back(gs);
  auto time2 = std::chrono::high_resolution_clock::now();
  auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
  Logging::Info("Found {} segments after Pruning in {}ms", new_segments.size(), total_time);

  return new_segments;
}

std::vector<GenomicSegment> PruneSegmentsByClusteringPerSeedSegment(
    const std::vector<GenomicSegment>& seed_segments,
    const std::vector<GenomicSegment>& segments_to_prune,
    const bool joint_clustering,
    std::optional<f64> h) {
  std::vector<GenomicSegment> res;
  IntervalTrees seg_trees(segments_to_prune);
  for (const auto& seed_seg : seed_segments) {
    std::vector<size_t> seg_idxs = seg_trees.LookUp(seed_seg.contig, seed_seg.start, seed_seg.end);
    if (!seg_idxs.empty()) {
      std::vector<GenomicSegment> segs;
      segs.reserve(seg_idxs.size());
      for (const auto& idx : seg_idxs) {
        segs.emplace_back(segments_to_prune[idx]);
      }
      std::vector<GenomicSegment> pruned_segs = PruneSegmentsByClustering(segs, joint_clustering, h);
      res.insert(res.end(), pruned_segs.begin(), pruned_segs.end());
    } else {
      Logging::Debug(
          "No segments found for seed segment: {}:{}-{}. Skipping pruning for this one and removing this seed segment "
          "from the output...",
          seed_seg.contig,
          seed_seg.start,
          seed_seg.end);
    }
  }
  return res;
}
}  // namespace xoos::cnc::segmentation
