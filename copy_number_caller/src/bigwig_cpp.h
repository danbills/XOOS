#pragma once

#include <xoos/types/float.h>
extern "C" {
#include <libbigwig/bigWig.h>
}

#include <memory>
#include <string>

#include <xoos/types/fs.h>

namespace xoos::cnc {
using BigWigPointer = std::unique_ptr<bigWigFile_t, decltype(&bwClose)>;
using BwStatsPointer = std::unique_ptr<f64[], decltype(&std::free)>;
using BwOverlappingIntervalsPtr = std::unique_ptr<bwOverlappingIntervals_t, decltype(&bwDestroyOverlappingIntervals)>;

class BigWig {
 public:
  BigWig() = delete;

  ~BigWig();

  explicit BigWig(const fs::path& fname);

  f64 GetMean(const std::string& contig, size_t start, size_t end);
  BwOverlappingIntervalsPtr GetInterval(const std::string& contig, size_t start, size_t end);

 private:
  BigWigPointer _fp;
};

}  // namespace xoos::cnc
