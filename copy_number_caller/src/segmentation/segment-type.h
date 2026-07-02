#pragma once

namespace xoos::cnc::segmentation {
enum class SegmentType {
  kUnknown,
  kSeed,
  kLogROnly,
  kBaf,
  kGermlineLikelihood,
  kSomaticWithBafLikelihood,
  kSomaticNoBafLikelihood
};

}  // namespace xoos::cnc::segmentation
