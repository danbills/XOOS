#pragma once

#include <xoos/types/float.h>
#include <xoos/types/fs.h>

namespace xoos::cnc {

// Maximum LogR value for a segment to be considered
const f64 kSegMaxLogR = 0.9;
const size_t kSegMinNumLogRs = 11;
const size_t kSegMinNumSnps = 11;

struct PurityPloidySearchOptions {
  f64 seg_max_logr = kSegMaxLogR;
  size_t seg_min_num_logrs = kSegMinNumLogRs;
  size_t seg_min_num_snps = kSegMinNumSnps;
};
}  // namespace xoos::cnc
