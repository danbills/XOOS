#include "lut-bundle-simplex.h"

#include <utility>

#include "adapter-design/adapter-design-bundle.h"
#include "lut-bundle/lut-bundle.h"

namespace xoos::demux {

LutBundleSimplex::LutBundleSimplex(SeqMatcher fixed_1, SeqMatcher sid_2, SeqMatcher fixed_3, SeqMatcher fixed_5,
                                   SeqMatcher sid_6, SeqMatcher fixed_7)
    : fixed_1_matcher{std::move(fixed_1)},
      sid_2_matcher{std::move(sid_2)},
      fixed_3_matcher{std::move(fixed_3)},
      fixed_5_matcher{std::move(fixed_5)},
      sid_6_matcher{std::move(sid_6)},
      fixed_7_matcher{std::move(fixed_7)} {}

template <>
LutBundleSimplex CreateLutBundle<LutBundleSimplex>(const AdapterDesignBundle& designs) {
  using enum BarcodeType;
  return LutBundleSimplex{*designs.GetMatcher5p(kAnchor), *designs.GetMatcher5p(kSid), *designs.GetMatcher5p(kStem),
                          *designs.GetMatcher3p(kStem),   *designs.GetMatcher3p(kSid), *designs.GetMatcher3p(kAnchor)};
}

}  // namespace xoos::demux
