#include "demux-and-trim-pipeline.h"

#include <xoos/log/logging.h>

#include "metrics/metric-filenames.h"
#include "task/flow-manager.h"

namespace xoos::demux {

namespace fs = std::filesystem;

/**
 * Check metrics validation and delete any existing metrics files if overwrite is enabled. If overwrite is not enabled
 * and any metrics files already exist, an error will be thrown.
 * Because metrics are written in the FlowManager destructor, we need to check for existing metrics files and handle
 * them before the FlowManager is created. This function is called in the DemuxAndTrimPipeline function before the
 * FlowManager is initialized.
 * @param param
 */
void HandleMetricsFiles(const DemuxAndTrimParam& param) {
  const auto metrics_dir = param.out_dir / kMetricsDirectory;
  if (!fs::exists(metrics_dir)) {
    fs::create_directories(metrics_dir);
    return;
  }
  if (!fs::is_directory(metrics_dir)) {
    throw error::Error(fmt::format(
        "Metrics path {} already exists but is not a directory. Please delete or rename this path, or choose a "
        "different output directory.",
        metrics_dir.string()));
  }
  // if overwriting, then by default each metrics overwrites the existing file within the code
  if (param.overwrite) {
    return;
  }
  const auto metrics_files = std::vector{
      // duplex metrics file names used in demux-and-trim pipeline
      metrics_dir / kRunMetricsFile, metrics_dir / kSampleMetricsFile, metrics_dir / kSampleAssignmentMetrics,

      // read length metrics file names
      metrics_dir / kUntrimmedFullReadLenDist, metrics_dir / kUntrimmedPartialReadLenDist,
      metrics_dir / kTrimmedFullReadLenDist, metrics_dir / kTrimmedPartialReadLenDist,

      // metrics file names
      metrics_dir / kPassingReadLengthDistr, metrics_dir / kFullDuplexReadLengthDistr,
      metrics_dir / kPartialDuplexReadLengthDistr, metrics_dir / kEndadapterPositionDistr,
      metrics_dir / kTotalReadLengthDistr, metrics_dir / kUnassignedReadLengthDistr,
      metrics_dir / kNoHairpinReadLengthDistr};
  for (const auto& file : metrics_files) {
    if (fs::exists(file)) {
      throw error::Error(
          fmt::format("Metrics file {} already exists. Please delete it or specify the --overwrite flag to allow "
                      "overwriting existing metrics files.",
                      file.string()));
    }
  }
}

void DemuxAndTrimPipeline(const DemuxAndTrimParam& param) {
  const auto& design = LoadAdapterDesign(param.adapter_design_bundle, param.adapter_design_name);
  fs::create_directories(param.out_dir);
  HandleMetricsFiles(param);
  const FlowManager flow_manager(param, design);
}

}  // namespace xoos::demux
