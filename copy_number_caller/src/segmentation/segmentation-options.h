#pragma once

#include <xoos/types/float.h>
#include <xoos/types/fs.h>

#include "segmentation/max-t.h"
#include "segmentation/segmentation-mode.h"

namespace xoos::cnc {
using segmentation::CbsMaxTMethod;
const f64 kSegmentationDefaultMaxP = 0.01;
const size_t kSegmentationDefaultMinObsPerSegment = 2;
const size_t kSegmentationDefaultNPermutations = 10000;
const size_t kSegmentationDefaultMaxNumSegments = 500;
const CbsMaxTMethod kSegmentationDefaultCbsMethod = CbsMaxTMethod::kBruteForce;
const f64 kSegmentationDefaultPruningClusteringParameter = 0.25;
const f64 kSegmentationDefaultMinTForAutomaticSegmentation = 7;
const bool kSegmentationDefaultUndoLogrSegments = false;
const f64 kSegmentationDefaultUndoLogrSegmentsSdFactor = 3.0;
const f64 kSegmentationDefaultIncrementUndoLogrSegmentsSdFactor = 1.0;
const bool kSegmentationDefaultUndoDhSegments = false;
const f64 kSegmentationDefaultUndoDhSegmentsSdFactor = 3.0;
const f64 kSegmentationDefaultIncrementUndoDhSegmentsSdFactor = 1.0;

struct SegmentationOptions {
  f64 max_p = kSegmentationDefaultMaxP;
  size_t min_obs_per_segment = kSegmentationDefaultMinObsPerSegment;
  size_t n_permutations = kSegmentationDefaultNPermutations;
  bool disable_hierarchical_pruning = false;
  bool disable_merging = false;
  bool no_single_obvs_subsegments = true;
  bool truncate_outliers = false;
  CbsMaxTMethod cbs_method = kSegmentationDefaultCbsMethod;
  segmentation::SegmentationMode segmentation_mode{segmentation::SegmentationMode::kUnknown};
  f64 pruning_clustering_parameter = kSegmentationDefaultPruningClusteringParameter;
  f64 min_t_for_automatic_segmentation = kSegmentationDefaultMinTForAutomaticSegmentation;
  bool undo_logr_segments = kSegmentationDefaultUndoLogrSegments;
  f64 undo_logr_segments_sd_factor = kSegmentationDefaultUndoLogrSegmentsSdFactor;
  f64 increment_undo_logr_segments_sd_factor = kSegmentationDefaultIncrementUndoLogrSegmentsSdFactor;
  bool undo_dh_segments = kSegmentationDefaultUndoDhSegments;
  f64 undo_dh_segments_sd_factor = kSegmentationDefaultUndoDhSegmentsSdFactor;
  f64 increment_undo_dh_segments_sd_factor = kSegmentationDefaultIncrementUndoDhSegmentsSdFactor;
  size_t max_num_segments = kSegmentationDefaultMaxNumSegments;
};
}  // namespace xoos::cnc
