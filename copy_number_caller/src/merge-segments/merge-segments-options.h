#pragma once

#include <xoos/types/float.h>
#include <xoos/types/fs.h>

namespace xoos::cnc {
const f64 kMergeDefaultMinMapq = 30;
const bool kMergeSegmentsDefaultUseMapqsObservations = false;
const bool kMergeSegmentsDefaultRecalculatePerSegmentData = false;

struct MergeSegmentsOptions {
  f64 min_mapq_threshold{};
  bool use_mapqs_observations = kMergeSegmentsDefaultUseMapqsObservations;
  bool recalculate_per_segment_data = kMergeSegmentsDefaultRecalculatePerSegmentData;
};
}  // namespace xoos::cnc
