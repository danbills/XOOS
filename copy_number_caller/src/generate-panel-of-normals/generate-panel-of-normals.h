#pragma once

#include <copy-number-caller/copy-number-caller-options.h>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "panel-of-normals.h"

namespace xoos::cnc {

struct GeneratePanelOfNormalsOut {
  PanelOfNormals on_target_reference_panel;
  PanelOfNormals off_target_reference_panel;
};

s32 GeneratePanelOfNormalsMain(const CopyNumberCallerOptions& options);
GeneratePanelOfNormalsOut GeneratePanelOfNormals(const CopyNumberCallerOptions& options);
}  // namespace xoos::cnc
