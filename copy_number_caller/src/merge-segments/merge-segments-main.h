#pragma once
#include <xoos/types/int.h>

#include "copy-number-caller/copy-number-caller-options.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {

using segmentation::GenomicSegment;
void MergeSegmentsMain(const CopyNumberCallerOptions& options);
}  // namespace xoos::cnc
