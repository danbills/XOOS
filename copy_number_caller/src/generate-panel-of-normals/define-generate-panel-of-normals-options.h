#pragma once
#include <copy-number-caller/copy-number-caller-options.h>

#include <xoos/cli/cli.h>
#include <xoos/types/fs.h>

namespace xoos::cnc {
void DefineOptionsGeneratePanelOfNormals(cli::AppPtr app, CopyNumberCallerOptions& options);
}  // namespace xoos::cnc
