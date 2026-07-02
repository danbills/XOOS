#pragma once
#include <xoos/types/fs.h>

#include "baits.h"
#include "copy-number-caller/copy-number-caller-options.h"

namespace xoos::cnc {

BaitRecords LoadOrGenerateWholeGenomeIntervals(const CopyNumberCallerOptions& options,
                                               const fs::path& augmented_baits_out);

}  // namespace xoos::cnc
