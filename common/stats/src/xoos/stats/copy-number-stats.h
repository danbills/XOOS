#pragma once

#include <xoos/types/float.h>

namespace xoos::stats {
f64 ExpectedLogR(f64 purity, f64 tumor_ploidy, f64 total_copy_number);
f64 ExpectedLogR(f64 purity, f64 tumor_ploidy, f64 total_copy_number, bool expect_haploid);
f64 ExpectedTotalCopyNumber(f64 purity, f64 tumor_ploidy, f64 mean_logr, bool expect_haploid);
f64 ExpectedAF(f64 purity, f64 copy_number, f64 multiplicity);
f64 ExpectedAF(f64 purity, f64 copy_number, f64 multiplicity, bool expect_haploid);
}  // namespace xoos::stats
