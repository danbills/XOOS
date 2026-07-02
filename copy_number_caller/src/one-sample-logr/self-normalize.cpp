#include "one-sample-logr/self-normalize.h"

#include <armadillo>
#include <cmath>

#include <xoos/log/logging.h>

#include "utility/utility-functions.h"

namespace xoos::cnc {

// purity=0.99, ploidy=2, expected logr formula: math.log2(2*(1-0.99)+(0.99*0)) / (2*(1-0.99) + (0.99*2))
const f64 kExpectedLogRForGermlineAtTCN0 = -6.643856;

/**
 * @brief self-normalize coverage counts by dividing by the median count on autosomes, and log2-transforming
 * @param coverage coverage records with count field populated
 * @return Observations object with logR values
 */
Observations SelfNormalizeCounts(const CoverageRecords& coverage) {
  // get the median count for the autosomes. TODO definitely can save memory by converting into a loop
  CoverageRecords autosome_coverage = coverage.GetAutosomes();
  f64 median_count = arma::median(autosome_coverage.count);
  if (median_count <= 0) {
    Logging::Error("median count for corrected coverage is 0!");
    throw std::runtime_error("no coverage!!");
  }
  Observations ret(coverage.count.size());
  for (size_t i = 0; i < coverage.count.size(); ++i) {
    auto [contig, start, end] = ParseRegionString(coverage.region[i]);
    ret.regions[i] = coverage.region[i];
    ret.contigs[i] = contig;
    ret.starts[i] = start;
    ret.ends[i] = end;
    if (coverage.count[i] > 0) {
      ret.obvs[i] = std::log2(coverage.count[i]) - std::log2(median_count);
      //  targets with 1 count can have logR value lower than the minimum we should expect; so just truncate these
      if (ret.obvs[i] < kExpectedLogRForGermlineAtTCN0) {
        ret.obvs[i] = kExpectedLogRForGermlineAtTCN0;
      }
    } else {
      ret.obvs[i] = kExpectedLogRForGermlineAtTCN0;
    }
  }
  return ret;
}

}  // namespace xoos::cnc
