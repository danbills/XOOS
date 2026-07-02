#include "demux-and-trim-duplex.h"

namespace xoos::demux {

struct DuplexMetrics;

static TrimDuplex CreateTrim(const LutBundleDuplex& lut_bundle) {
  return TrimDuplex{lut_bundle.Runway5pMatcher(), lut_bundle.StartSequenceMatcher()};
}

DemuxAndTrimDuplex::DemuxAndTrimDuplex(const LutBundleDuplex& lut_bundle) : _trim{CreateTrim(lut_bundle)} {}

}  // namespace xoos::demux
