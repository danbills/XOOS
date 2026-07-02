#pragma once
#include <armadillo>
#include <vector>

#include "segmentation/genomic-segments.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc::segmentation {
std::vector<Segment> CircularBinarySegmentationWrapper(const arma::vec& obvs, const SegmentationOptions& options);
}
