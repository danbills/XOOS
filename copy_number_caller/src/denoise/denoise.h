#pragma once
#include <string>
#include <tuple>
#include <vector>

#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "copy-number-caller/copy-number-caller-options.h"
#include "coverage.h"
#include "observations.h"
#include "panel-of-normals.h"
#include "sex.h"

namespace xoos::cnc {
using ObservationsWOnTarget = std::tuple<Observations, std::vector<bool>>;

struct DenoiseOut {
  Observations logrs;
  Observations denoised_logrs;
  std::vector<bool> on_target_status;
  Sex sex;
};

ObservationsWOnTarget CombineOnOffTargetLogRs(const std::vector<std::string>& regions,
                                              const Observations& on_target_logrs,
                                              const Observations& off_target_logrs);
void WriteSex(char sex, std::ostream& os);
Observations FilterLogrs(const std::vector<std::string>& intervals_to_keep, char sex, const Observations& logrs);
DenoiseOut Denoise(const CoverageRecords& tumor_cov,
                   const PanelOfNormals& on_target_ref_pool,
                   const PanelOfNormals& off_target_ref_pool,
                   bool no_filter,
                   size_t min_target_length,
                   f64 min_panel_median_cov,
                   f64 min_panel_median_and_tumor_cov,
                   f64 min_off_target_filter_frac);
void DenoiseMain(const CopyNumberCallerOptions& options);

}  // namespace xoos::cnc
