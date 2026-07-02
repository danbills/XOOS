#pragma once

#include <htslib/sam.h>

#include <xoos/io/alignment-reader.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "baits.h"
#include "calculate-coverage/calculate-coverage-options.h"
#include "copy-number-caller/copy-number-caller-options.h"
#include "coverage.h"

namespace xoos::cnc {

/**
 * Non-owning view passed to pileup callbacks. All pointed-to resources are
 * owned by the surrounding AlignmentReader / HtsItrPtr in CalculateCoverageRegion.
 */
struct Plpconf {
  htsFile* fp{nullptr};
  hts_itr_t* itr{nullptr};
  u32 exclude_flags{0};
  s32 total_reads{0};
  s32 sum_mapping_quality{0};
  bool marked_dup{false};
  bool had_read_error{false};
  bool has_rescued_secondaries{false};
};

struct CalculateCoverageOut {
  CoverageRecords cov;
};

// this will handle all the IO logic and load data structures to be passed into CalculateCoverage.
void CalculateCoverageMain(const CopyNumberCallerOptions& options);
CoverageRecords CalculateCoverage(const BaitRecords& baits, const CalculateCoverageOptions& options);
void CalculateCoverageRegion(const io::AlignmentReader& reader,
                             CoverageRecords& res,
                             u32 exclude_flags,
                             bool ignore_dn,
                             s32 start_idx,
                             s32 step);
}  // namespace xoos::cnc
