#include "segmentation/cbs-wrapper.h"

#include "segmentation/cbs-for-somatic.h"
#include "segmentation/cbs.h"

namespace xoos::cnc::segmentation {
std::vector<Segment> CircularBinarySegmentationWrapper(const arma::vec& obvs, const SegmentationOptions& options) {
  if (options.segmentation_mode == SegmentationMode::kGermline) {
    return CircularBinarySegmentation(obvs, options);
  } else if (options.segmentation_mode == SegmentationMode::kSomatic) {
    return somatic::CircularBinarySegmentation(
        obvs, options.max_p, options.n_permutations, options.min_obs_per_segment, options.no_single_obvs_subsegments);
  } else {
    throw std::runtime_error("invalid segmentation mode");
  }
}
}  // namespace xoos::cnc::segmentation
