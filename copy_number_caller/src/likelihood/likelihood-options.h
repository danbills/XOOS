#pragma once

#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "likelihood/likelihood-mode.h"
#include "likelihood/likelihood-model.h"
#include "segmentation/segment-type.h"

namespace xoos::cnc {
using segmentation::SegmentType;
const s32 kLikelihoodDefaultSomaticMAPQCutoffForCalls = -1;
const s32 kLikelihoodDefaultGermlineMAPQCutoffForCalls = 30;
const size_t kLikelihoodCnvLengthFlagMinSize = 2000;

struct LikelihoodOptions {
  LikelihoodModel likelihood_model = LikelihoodModel::kSerialSummarized;
  SegmentType input_segments_type{SegmentType::kUnknown};
  SegmentType output_segments_type{SegmentType::kUnknown};
  size_t cnv_length_flag_min_size = kLikelihoodCnvLengthFlagMinSize;
  s32 mapq_cutoff_for_calls{};
  bool mapq_cutoff_for_calls_is_user_set{false};
  LikelihoodMode mode{LikelihoodMode::kUnknown};
};

}  // namespace xoos::cnc
