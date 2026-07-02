#pragma once
#include <xoos/cli/cli.h>

#include <memory>

#include "core/demux-and-trim-pipeline.h"

namespace xoos::demux {

using DemuxAndTrimParamPtr = std::shared_ptr<DemuxAndTrimParam>;

// Define the command line options and wire them up to the parameter object.
[[maybe_unused]] void DefineOptions(cli::AppPtr app, const DemuxAndTrimParamPtr& param);

cli::PreCallback<DemuxAndTrimParam> CreatePreCallback();

/// Expands the input file list by recursively finding all sequence files in the directories specified in the input file
/// list. If an input file is not a directory, it is included in the expanded input file list as is.
std::vector<fs::path> ExpandInputFileList(const std::vector<fs::path>& input_files);
}  // namespace xoos::demux
