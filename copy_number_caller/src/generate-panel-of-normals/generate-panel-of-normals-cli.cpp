#include "generate-panel-of-normals/generate-panel-of-normals-cli.h"

#include <xoos/cli/cli.h>
#include <xoos/log/logging.h>

#include "copy-number-caller/copy-number-caller-options.h"
#include "generate-panel-of-normals/define-generate-panel-of-normals-options.h"
#include "generate-panel-of-normals/generate-panel-of-normals.h"

namespace xoos::cnc {

s32 GeneratePanelOfNormalsCliMain(s32 argc, char** argv) {
  auto app = cli::SetupDefaultCli(PROGRAM_NAME, VERSION);
  CopyNumberCallerOptions options;
  DefineOptionsGeneratePanelOfNormals(app.get(), options);
  app->require_subcommand(0);  /// require exactly 0 subcommand
  // this value will be returned by the main function and updated by whatever submodule is called
  s32 ret_val = EXIT_SUCCESS;
  app->callback([&options, &ret_val]() { ret_val = GeneratePanelOfNormalsMain(options); });
  CLI11_PARSE(*app, argc, argv);
  return ret_val;
}
}  // namespace xoos::cnc
