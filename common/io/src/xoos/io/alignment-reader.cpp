#include "xoos/io/alignment-reader.h"

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-util.h>

namespace xoos::io {

AlignmentReader OpenAlignmentReader(const fs::path& location) {
  auto bam = io::HtsOpen(location, "rb");
  if (bam->format.format != htsExactFormat::sam && bam->format.format != htsExactFormat::bam) {
    throw error::Error(
        "Unsupported alignment file format when opening '{}': detected HTS format code {} (supported formats: SAM, "
        "BAM)",
        location.string(),
        std::to_string(static_cast<int>(bam->format.format)));
  }

  auto idx = io::HtsIdxLoad(location, HTS_FMT_BAI);
  auto hdr = io::SamHdrRead(bam.get());
  return AlignmentReader{std::move(bam), std::move(hdr), std::move(idx)};
}

vec<AlignmentReader> OpenAlignmentReaders(const fs::path& bam_filename, const u32 reader_count) {
  vec<AlignmentReader> alignment_readers;
  alignment_readers.reserve(reader_count);
  for (u32 i = 0; i < reader_count; ++i) {
    alignment_readers.emplace_back(OpenAlignmentReader(bam_filename));
  }
  return alignment_readers;
}

}  // namespace xoos::io
