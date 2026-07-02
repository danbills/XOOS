#pragma once

namespace xoos::svc {
// Enum for including or excluding `duplex_lowbq` in calculated BAM features `duplex_dp` and `duplex_af`
enum class DuplexLowbqMode {
  kInclude,
  kExclude
};
}  // namespace xoos::svc
