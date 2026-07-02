#pragma once
#include <optional>

#include "copy-number-caller/common/coverage-check.h"
#include "likelihood/likelihood.h"
#include "likelihood/total-copy-number-prior.h"
#include "observations.h"
#include "sex.h"
#include "vcf-purity-source.h"

namespace xoos::cnc {
struct CopyNumberCallerRet {
  Observations ref_depths;
  Observations alt_depths;
  std::vector<GenomicSegment> segments{};
  Sex sex{};
  std::optional<f64> tumor_purity;
  std::optional<f64> tumor_ploidy;
  std::optional<VcfPuritySource> vcf_purity_source;
  std::optional<TotalCopyNumberPrior> total_copy_number_prior;
  std::optional<CoverageCheckResult> coverage_check_result;
};

}  // namespace xoos::cnc
