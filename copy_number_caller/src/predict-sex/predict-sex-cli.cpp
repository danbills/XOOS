#include "predict-sex/predict-sex-cli.h"

#include <predict-sex/define-predict-sex-options.h>

#include <CLI/CLI.hpp>

#include <xoos/cli/cli.h>
#include <xoos/log/logging.h>

#include "predict-sex/predict-sex-options.h"
#include "predict-sex/predict-sex.h"

namespace xoos::cnc {

s32 PredictSexCliMain(s32 argc, char** argv) {
  // this object handles all the logic for parsing command line arguments
  const auto app = cli::SetupDefaultCli(PROGRAM_NAME, VERSION);
  // this will store all our parameters
  PredictSexOptions options;
  // command line argument definitions are in this function, which also populate the above object
  DefineOptionsPredictSex(app.get(), options);
  app->require_subcommand(0, 1);  /// require exactly 0 or 1 subcommands
  // this value will be returned by the main function and updated by whatever submodule is called
  s32 ret_val = EXIT_SUCCESS;
  app->callback([&ret_val, &options]() { ret_val = PredictSexMain(options); });
  // this is a macro provided by CLI11 to parse the arguments. It raises an exception if it encounters an error during
  // parsing.
  CLI11_PARSE(*app, argc, argv);
  return ret_val;
}
}  // namespace xoos::cnc
