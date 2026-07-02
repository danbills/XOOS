#include "bigwig_cpp.h"

#include <xoos/error/error.h>

namespace xoos::cnc {
BigWig::BigWig(const fs::path& fname) : _fp(BigWigPointer{bwOpen(fname.c_str(), nullptr, "r"), bwClose}) {
  // not calling bwInit because it has to occur before assigning _fp, and we cannot initialize _fp in here since it's a
  // unique pointer. This should be fine because bwInit is only really needed for remote BigWig connections, which I
  // don't think will happen here
  if (!_fp) {
    throw std::runtime_error(
        "The provided BigWig (.bw) file is invalid due to corruption, truncation, "
        "or unsupported format. Please check the file's integrity and format.");
  }
}

BigWig::~BigWig() {
  bwCleanup();
}

/**
 * @brief get mean value over the entire interval
 * @param contig
 * @param start
 * @param end
 * @return
 */
f64 BigWig::GetMean(const std::string& contig, size_t start, size_t end) {
  auto dbl_arr = BwStatsPointer(bwStats(_fp.get(), contig.c_str(), start, end, 1, mean), std::free);
  if (!dbl_arr) {
    throw error::Error("Could not retrieve region {}:{}-{} from BigWig", contig, start, end);
  }
  // dbl-arr should only have one element. this should copy the f64 over to the return value and then destroys the
  // ptr
  return dbl_arr.get()[0];
}

BwOverlappingIntervalsPtr BigWig::GetInterval(const std::string& contig, size_t start, size_t end) {
  return {bwGetOverlappingIntervals(_fp.get(), contig.data(), start, end), bwDestroyOverlappingIntervals};
}

}  // namespace xoos::cnc
