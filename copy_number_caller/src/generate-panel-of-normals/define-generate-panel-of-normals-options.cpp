#include "generate-panel-of-normals/define-generate-panel-of-normals-options.h"

#include <copy-number-caller/copy-number-caller-options.h>

#include <xoos/cli/cli.h>
#include <xoos/cli/validators/file-permission-validator.h>

#include "copy-number-caller/define-copy-number-caller-options.h"
#include "misc/define-thread-options.h"

namespace xoos::cnc {
void DefineOptionsGeneratePanelOfNormals(cli::AppPtr app, CopyNumberCallerOptions& options) {
  DefineAugmentBaitsInputOptions(app, options, false);
  DefineIntervalOptions(app,
                        options.augment_baits_options,
                        false,
                        kAugmentBaitsDefaultWholeGenomeIntervalSize,
                        kAugmentBaitsDefaultWholeGenomeMinMappability);
  DefineCalculateCoverageAlgorithmOptions(app, options.calculate_coverage_options, false);
  DefineGCCorrectAlgorithmOptions(app, options.gc_correct_options, false);
  DefineThreadOptions(app, options.threads);
  // cli::AddThreadCountOption(app, opt::kThreads, options.threads)->group(opt::kPerformanceOptionsGroup);
  app->add_option("--panel-of-normals-bams",
                  options.panel_of_normals_lists,
                  "text file containing coverage file locations for panel of normals bams")
      ->check(cli::FileListReadableValidator())
      ->required();
  app->add_option("--output-dir",
                  options.output_dir,
                  "Directory where panel-of-normals outputs (including the HDF5 file) will be written");
  app->add_option("--panel-of-normals-hdf5",
                  options.panel_of_normals_hdf5_fname,
                  "Filename for the generated panel-of-normals HDF5 file (relative to --output-dir)");
}
}  // namespace xoos::cnc
