#include "one-sample-logr/one-sample-logr.h"

#include <xoos/util/file-functions.h>

#include "copy-number-caller/copy-number-caller-options.h"
#include "coverage.h"
#include "io/column-names.h"
#include "io/copy-number-caller-default-filenames.h"
#include "observations.h"
#include "one-sample-logr/self-normalize.h"

namespace xoos::cnc {
void OneSampleLogRMain(const CopyNumberCallerOptions& options) {
  const auto logrs_out = options.output_dir / kDefaultLogRsOutput;
  file::CheckFilePermissionsAndOutputPathExistence({}, {logrs_out});
  std::ifstream cov_ifs(options.normal_coverage_fname.value());
  CoverageRecords cov(cov_ifs);
  Observations obvs(SelfNormalizeCounts(cov));
  std::ofstream ofs(logrs_out);
  obvs.Write(ofs, kColumnLogRatio, false, options.command_line_info);
}
}  // namespace xoos::cnc
