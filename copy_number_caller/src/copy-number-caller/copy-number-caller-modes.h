#pragma once

namespace xoos::cnc {
enum class CopyNumberCallerModes {
  kSomaticTumorNormalWGS,
  kSomaticTumorTargetedEnrichment,
  kGermlineNormalWGS,
  kUnknown
};
}  // namespace xoos::cnc
