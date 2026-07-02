#include "lut-bundle-simplex-10x.h"

#include <utility>

#include "adapter-design/adapter-design-bundle.h"
#include "lut-bundle/lut-bundle.h"

namespace xoos::demux {

LutBundleSimplex10x::LutBundleSimplex10x(SeqMatcher fixed_1,  // NOSONAR(S5419) -- mirrors LutBundleSimplex convention
                                         SeqMatcher sid_2, SeqMatcher fixed_3)
    : fixed_1_matcher{std::move(fixed_1)}, sid_2_matcher{std::move(sid_2)}, fixed_3_matcher{std::move(fixed_3)} {}

template <>
LutBundleSimplex10x CreateLutBundle<LutBundleSimplex10x>(
    const AdapterDesignBundle& designs) {  // NOSONAR(S5536) -- called from flow-manager
  using enum BarcodeType;
  return LutBundleSimplex10x{*designs.GetMatcher5p(kAnchor), *designs.GetMatcher5p(kSid), *designs.GetMatcher5p(kStem)};
}

}  // namespace xoos::demux
