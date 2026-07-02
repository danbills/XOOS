#pragma once

#include <xoos/cli/cli.h>

#include "rescue-secondary/rescue-secondary-options.h"

namespace xoos::read_collapser::rescue_secondary {

void SetRescueSecondaryCommandLineInfo(cli::ConstAppPtr app, const RescueSecondaryOptionsPtr& options);

void DefineRescueSecondaryOptions(cli::AppPtr app, const RescueSecondaryOptionsPtr& options);

}  // namespace xoos::read_collapser::rescue_secondary
