#include "generate-panel-of-normals/generate-panel-of-normals.h"

#include <copy-number-caller/copy-number-caller-options.h>

#include <taskflow/algorithm/for_each.hpp>

#include <xoos/log/logging.h>
#include <xoos/util/file-functions.h>

#include "augment-baits/augment-baits.h"
#include "copy-number-caller/common/bam-to-corrected-coverage.h"
#include "io/copy-number-caller-default-filenames.h"
#include "utility/utility-functions.h"

namespace xoos::cnc {
s32 GeneratePanelOfNormalsMain(const CopyNumberCallerOptions& options) {
  GeneratePanelOfNormalsOut out = GeneratePanelOfNormals(options);
  Logging::Info("writing panel of normals to {}", options.panel_of_normals_hdf5_fname.string());
  // out.on_target_reference_panel.SerializeToHDF5(options.panel_of_normals_hdf5_fname, false);
  // out.off_target_reference_panel.SerializeToHDF5(options.panel_of_normals_hdf5_fname, true);
  return EXIT_SUCCESS;
}

GeneratePanelOfNormalsOut GeneratePanelOfNormals(const CopyNumberCallerOptions& options) {
  BaitRecords baits = GenerateAndAugmentOnAndOffTargetBaits(options);
  const auto augmented_baits_out = options.output_dir / kDefaultAugmentedBaitsOutput;
  file::CheckFilePermissionsAndOutputPathExistence({}, {augmented_baits_out});
  std::ofstream ofs(augmented_baits_out);
  Logging::Info("Writing baits file to {}", augmented_baits_out.string());
  baits.Write(ofs, options.command_line_info);
  tf::Taskflow taskflow;
  // calculate coverage
  std::vector<fs::path> panel_of_normals_fnames(PathsFromFile(options.panel_of_normals_lists));
  std::vector<BamToCorrectedCoverageTaskFlowGraph> normal_corrected_coverage_tf_graphs;
  normal_corrected_coverage_tf_graphs.reserve(panel_of_normals_fnames.size());
  std::vector<tf::Task> panel_of_normals_coverage_and_correction_tasks(panel_of_normals_fnames.size());
  for (size_t i = 0; i < panel_of_normals_fnames.size(); ++i) {
    // We have to use emplace_back instead of assigning since the default empty constructor is deleted, and I think the
    // move constructor is also deleted as a result
    normal_corrected_coverage_tf_graphs.emplace_back(panel_of_normals_fnames[i],
                                                     baits,
                                                     std::nullopt,
                                                     std::nullopt,
                                                     options.calculate_coverage_options.exclude_flags,
                                                     options.calculate_coverage_options.ignore_DN,
                                                     options.gc_correct_options.first_span,
                                                     options.command_line_info);
    panel_of_normals_coverage_and_correction_tasks[i] = taskflow.composed_of(normal_corrected_coverage_tf_graphs[i]);
  }
  GeneratePanelOfNormalsOut res;
  tf::Task create_panel_task = taskflow.emplace([&normal_corrected_coverage_tf_graphs, &res]() {
    // collect corrected normal coverages
    std::vector<CoverageRecords> normal_covs;
    normal_covs.reserve(normal_corrected_coverage_tf_graphs.size());
    for (auto& tf_graph : normal_corrected_coverage_tf_graphs) {
      normal_covs.emplace_back(tf_graph.GetResult());
    }
    if (std::any_of(
            normal_covs.begin(), normal_covs.end(), [](const CoverageRecords& cov) { return cov.count.empty(); })) {
      throw std::runtime_error("panel of normals sample has empty coverage!");
    }
    Logging::Info("building on-target panel of normals");
    res.on_target_reference_panel = PanelOfNormals(normal_covs, true);
    Logging::Info("building off-target panel of normals");
    res.off_target_reference_panel = PanelOfNormals(normal_covs, false);
  });
  // define taskflow graph
  for (auto& task : panel_of_normals_coverage_and_correction_tasks) {
    task.precede(create_panel_task);
  }
  Logging::Info("Using {} threads", options.threads);
  tf::Executor executor(options.threads);
  executor.run(taskflow).get();
  return res;
}
}  // namespace xoos::cnc
