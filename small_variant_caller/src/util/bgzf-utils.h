#pragma once

#include <filesystem>
#include <memory>

#include <htslib/bgzf.h>

#include <xoos/types/vec.h>

namespace xoos::svc {

namespace fs = std::filesystem;

/**
 * @brief RAII deleter for BGZF* handles.
 *
 * Calls bgzf_close on destruction, flushing any buffered data and writing the BGZF EOF marker.
 * Used as the custom deleter for BgzfPtr.
 */
struct BgzfCloser {
  void operator()(BGZF* const f) const {
    if (f != nullptr) {
      bgzf_close(f);
    }
  }
};

/// Owning smart pointer for a BGZF handle. Calls bgzf_close on destruction.
using BgzfPtr = std::unique_ptr<BGZF, BgzfCloser>;

/**
 * @brief Concatenate a header BGZF file and per-region BGZF data files into a single output.
 *
 * Performs block-level BGZF concatenation (no decompression/recompression). Per-region files
 * must contain only BGZF-compressed VCF record lines with no VCF header.
 *
 * Layout of the final output:
 *   [header BGZF blocks] [region_0 data blocks] ... [region_N data blocks] [EOF marker]
 *
 * @param output_path Destination file path for the concatenated VCF.
 * @param header_file BGZF file containing only the VCF header.
 * @param region_files Ordered list of per-region BGZF data files to append after the header.
 * @throws std::runtime_error if any file cannot be opened or a read/write error occurs.
 */
void ConcatenateBgzfFiles(const fs::path& output_path, const fs::path& header_file, const vec<fs::path>& region_files);

}  // namespace xoos::svc
