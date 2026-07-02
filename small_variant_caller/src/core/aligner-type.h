#pragma once

namespace xoos::svc {

// Aligner type for tumor-normal-wgs model selection.
// kCustom requires the user to explicitly provide --model and score thresholds.
enum class AlignerType {
  kBwa,
  kGiraffe,
  kCustom
};

}  // namespace xoos::svc
