#pragma once
#include <vector>

#include "segmentation/genomic-segments.h"
#include "segmentation/segment-type.h"
#include "segmentation/segments-header.h"

namespace xoos::cnc::segmentation {

std::vector<GenomicSegment> ReadSegments(const fs::path& fname, SegmentType segment_type);
SegmentsHeader ReadHeaderFromSegments(const fs::path& fname);
}  // namespace xoos::cnc::segmentation
