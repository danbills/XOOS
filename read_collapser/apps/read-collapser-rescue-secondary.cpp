#include <xoos/cli/cli.h>

#include "rescue-secondary/rescue-secondary-cli.h"
#include "rescue-secondary/rescue-secondary-options.h"
#include "rescue-secondary/rescue-secondary.h"

int main(int argc, char** argv) {
  xoos::cli::StandardMainParam<xoos::read_collapser::rescue_secondary::RescueSecondaryOptions> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<xoos::read_collapser::rescue_secondary::RescueSecondaryOptions>(),
      .define_options = xoos::read_collapser::rescue_secondary::DefineRescueSecondaryOptions,
      .main = xoos::read_collapser::rescue_secondary::RescueSecondary,
      .pre_callback = xoos::read_collapser::rescue_secondary::SetRescueSecondaryCommandLineInfo,
  };
  return StandardMain(argc, argv, param);
}
