#pragma once

#include "rescue-secondary/rescue-secondary-options.h"

namespace xoos::read_collapser::rescue_secondary {

// Entry point for the rescue secondary alignment workflow.
// Dispatches to collated or non-collated processing based on options.
void RescueSecondary(const RescueSecondaryOptions& options);

}  // namespace xoos::read_collapser::rescue_secondary
