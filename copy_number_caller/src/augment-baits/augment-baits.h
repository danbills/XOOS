#pragma once
#include <xoos/types/int.h>

#include "baits.h"
#include "copy-number-caller/copy-number-caller-options.h"

namespace xoos::cnc {
void AugmentBaitsMain(const CopyNumberCallerOptions& options);
BaitRecords GenerateAndAugmentOnAndOffTargetBaits(const CopyNumberCallerOptions& options);
BaitRecords GenerateAndAugmentWholeGenomeIntervals(const CopyNumberCallerOptions& options);
}  // namespace xoos::cnc
