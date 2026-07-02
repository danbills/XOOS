#pragma once

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

namespace xoos::read_collapser::rescue_secondary {

constexpr std::string_view kDefaultMetricsFilename = "rescue_secondary_metrics.tsv";

struct RescueSecondaryMetrics {
  u64 num_records{};
  u64 num_primary_records{};
  u64 num_secondary_records{};
  u64 num_supplementary_records{};
  u64 num_rescued_secondary_records{};

  void WriteToTsv(const fs::path& output, const io::CommandLineInfo& command_line_info) const;
};

}  // namespace xoos::read_collapser::rescue_secondary
