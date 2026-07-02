#pragma once
#include <map>
#include <string>

namespace xoos::cnc::segmentation {
enum class SegmentationMode {
  kGermline,
  kSomatic,
  kUnknown
};
extern const std::map<std::string, SegmentationMode> kStringToSegmentationMode;
}  // namespace xoos::cnc::segmentation
