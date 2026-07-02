#pragma once
#include <vector>

#include <xoos/types/float.h>

#include "segmentation/genomic-segments.h"

namespace xoos::cnc {
using segmentation::GenomicSegment;
GenomicSegment MergeSegments(const std::vector<GenomicSegment>& segments_to_merge);
std::vector<GenomicSegment> MergeLowMapqSegments(const std::vector<GenomicSegment>& segments, f64 min_mapq_cutoff);
std::vector<GenomicSegment> MergeAdjacentEqualCopyNumberSegments(const std::vector<GenomicSegment>& segments);
}  // namespace xoos::cnc
