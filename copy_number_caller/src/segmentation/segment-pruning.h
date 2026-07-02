#pragma once

#include <optional>
#include <vector>

#include <xoos/types/float.h>

#include "observations.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc::segmentation {
std::vector<GenomicSegment> BreakLogrSegmentsByDHSegments(const std::vector<GenomicSegment>& logr_segments,
                                                          const Observations& logrs,
                                                          const std::vector<GenomicSegment>& dh_segments,
                                                          const Observations& dhs);
std::vector<GenomicSegment> PruneSegmentsByClustering(const std::vector<GenomicSegment>& logr_segments,
                                                      bool joint_clustering,
                                                      std::optional<f64> h);
std::vector<GenomicSegment> PruneSegmentsByClusteringPerSeedSegment(
    const std::vector<GenomicSegment>& seed_segments,
    const std::vector<GenomicSegment>& segments_to_prune,
    bool joint_clustering,
    std::optional<f64> h);
}  // namespace xoos::cnc::segmentation
