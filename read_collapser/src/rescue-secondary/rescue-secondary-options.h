#pragma once

#include <memory>
#include <optional>
#include <string>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

namespace xoos::read_collapser::rescue_secondary {

constexpr std::string kStdin = "-";

struct RescueSecondaryOptions {
  // Input options
  std::string bam{};
  bool collated{false};

  // Output options
  fs::path output_dir{};
  std::optional<std::string> prefix{};
  bool overwrite{false};

  // Rescue options
  f64 min_alignment_score_ratio{};

  // Performance options
  size_t threads{};

  // Program metadata
  io::CommandLineInfo command_line_info{};
};

using RescueSecondaryOptionsPtr = std::shared_ptr<RescueSecondaryOptions>;

}  // namespace xoos::read_collapser::rescue_secondary
