#include "gc-correct/gc-correct.h"

#include <fstream>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/gc_correct/gc-correct-taskflow-graph.h>
#include <xoos/log/logging.h>
#include <xoos/util/file-functions.h>

#include "baits.h"
#include "copy-number-caller/copy-number-caller-options.h"
#include "coverage.h"
#include "io/copy-number-caller-default-filenames.h"

namespace xoos::cnc {
void GCCorrectMain(const CopyNumberCallerOptions& options) {
  const auto normal_corrected_coverage_out = options.output_dir / kDefaultNormalCorrectedCoverageOutput;
  file::CheckFilePermissionsAndOutputPathExistence({}, {normal_corrected_coverage_out});
  tf::Taskflow taskflow;
  tf::Executor executor(options.threads);
  Logging::Info("Loading and GC correcting coverage for sample {}", options.normal_coverage_fname.value().string());
  std::ifstream cov_istream(options.normal_coverage_fname.value());
  std::ifstream bait_istream(options.augmented_baits_fname.value());
  CoverageRecords covs(cov_istream);
  BaitRecords baits = BaitRecords(bait_istream);
  if (!VerifyRegions(covs.region, baits.GetRegions())) {
    throw error::Error("Regions in coverage and baits files do not match");
  }
  gc_correct::GCCorrectTaskFlowGraph gc_correction_taskflow_graph(covs.region,
                                                                  covs.count,
                                                                  covs.total_coverage,
                                                                  baits.GetGCBias(),
                                                                  baits.GetMappability(),
                                                                  baits.GetOnTargetStatus(),
                                                                  options.gc_correct_options.first_span);
  taskflow.composed_of(gc_correction_taskflow_graph);  // NOLINT
  executor.run(taskflow).get();
  gc_correct::GCCorrectResults res = gc_correction_taskflow_graph.GetResult();
  CoverageRecords corrected_covs;
  corrected_covs.region = covs.region;
  corrected_covs.total_coverage = res.total_coverage;
  corrected_covs.average_coverage = res.average_coverage;
  corrected_covs.count = res.counts;
  corrected_covs.on_target = covs.on_target;
  Logging::Info("Printing tumor counts to {}", normal_corrected_coverage_out.string());
  std::ofstream ofs(normal_corrected_coverage_out);
  corrected_covs.Write(ofs, options.command_line_info);
}
}  // namespace xoos::cnc
