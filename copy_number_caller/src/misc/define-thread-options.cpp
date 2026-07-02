#include "misc/define-thread-options.h"

#include <CLI/CLI.hpp>

#include <xoos/cli/cli.h>
#include <xoos/cli/thread-count-option-util.h>

namespace xoos::cnc {

CLI::Option_group* DefineThreadOptions(cli::AppPtr app, size_t& threads) {
  auto* option_group = app->add_option_group("Thread");
  cli::AddThreadCountOption(option_group, "--threads", threads);
  return option_group;
}

}  // namespace xoos::cnc
