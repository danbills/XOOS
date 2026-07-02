#pragma once
#include <xoos/cli/cli.h>

#include "predict-sex/predict-sex-options.h"

namespace xoos::cnc {

void DefineOptionsPredictSex(cli::AppPtr app, PredictSexOptions& options);
}  // namespace xoos::cnc
