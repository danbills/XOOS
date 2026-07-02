#include "segmentation/segmentation-mode.h"

namespace xoos::cnc::segmentation {
const std::map<std::string, SegmentationMode> kStringToSegmentationMode{
    {"germline", SegmentationMode::kGermline},
    {"somatic", SegmentationMode::kSomatic},
};
}  // namespace xoos::cnc::segmentation
