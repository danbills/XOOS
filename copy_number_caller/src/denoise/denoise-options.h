#pragma once

#include <xoos/types/float.h>
#include <xoos/types/fs.h>

namespace xoos::cnc {
const size_t kDenoiseDefaultMinTargetLength = 5;
const f64 kDenoiseDefaultMinPanelMedianCov = 15;
const f64 kDenoiseDefaultMinPanelMedianAndTumorCov = 100;
const f64 kDenoiseDefaultMinOffTargetFilterFrac = 0.05;
const bool kDenoiseDefaultNoFilter = false;

struct DenoiseOptions {
  size_t min_target_length = 5;
  f64 min_panel_median_cov = 15;
  f64 min_panel_median_and_tumor_cov = 100;
  f64 min_off_target_filter_frac = 0.05;
  bool no_filter = false;
};
}  // namespace xoos::cnc
