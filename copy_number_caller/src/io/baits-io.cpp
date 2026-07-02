#include "io/baits-io.h"

#include <fstream>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include "augment-baits/augment-baits.h"

namespace xoos::cnc {

BaitRecords LoadOrGenerateWholeGenomeIntervals(const CopyNumberCallerOptions& options,
                                               const fs::path& augmented_baits_out) {
  if (options.augmented_baits_fname.has_value()) {
    std::ifstream ifs(options.augmented_baits_fname.value());
    return BaitRecords(ifs);
  } else {
    auto baits = GenerateAndAugmentWholeGenomeIntervals(options);
    Logging::Info("Writing augmented baits to: {}", augmented_baits_out.string());
    std::ofstream ofs(augmented_baits_out);
    if (!ofs.is_open()) {
      throw error::Error("Failed to open augmented baits output file: {}", augmented_baits_out.string());
    }
    baits.Write(ofs, options.command_line_info);
    return baits;
  }
}

}  // namespace xoos::cnc
