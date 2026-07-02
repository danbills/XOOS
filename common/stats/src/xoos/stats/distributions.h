#pragma once

#include <xoos/types/float.h>

namespace xoos::stats {
f64 BinomialPMF(f64 k, f64 n, f64 p);
f64 LogBinomialPMF(f64 k, f64 n, f64 p);
f64 NonStandardTCDF(f64 x, f64 df, f64 ncp);
f64 NonStandardTCDF(f64 x, f64 df, f64 mean, f64 sd);
f64 NonStandardTPDF(f64 x, f64 df, f64 ncp);
f64 NonStandardTPDF(f64 x, f64 df, f64 mean, f64 sd);
}  // namespace xoos::stats
