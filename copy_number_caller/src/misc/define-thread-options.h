#pragma once
#include <cstddef>

#include <CLI/CLI.hpp>

#include <xoos/cli/cli.h>

namespace xoos::cnc {
const size_t kThreadOptionsDefaultNThreads = 1;
CLI::Option_group* DefineThreadOptions(cli::AppPtr app, size_t& threads);
}  // namespace xoos::cnc
