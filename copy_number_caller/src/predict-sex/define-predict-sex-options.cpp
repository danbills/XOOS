#include "predict-sex/define-predict-sex-options.h"

#include <xoos/cli/validators/file-permission-validator.h>

namespace xoos::cnc {

void DefineOptionsPredictSex(cli::AppPtr app, PredictSexOptions& options) {
  app->add_option("--coverage-file", options.coverage_fname, "coverage file (pre or post GC correction)")
      ->check(cli::FileReadableValidator());
}

}  // namespace xoos::cnc
