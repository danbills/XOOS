#include "lut-bundle-ysu.h"

#include <xoos/error/error.h>

#include <utility>

#include "adapter-design/adapter-design-bundle.h"
// Clangd's include-cleaner doesn't detect usage from template specializations.
// This include is required for explicit specialization below.
#include "lut-bundle/lut-bundle.h"  // IWYU pragma: keep

namespace xoos::demux {

LutBundleYsu::LutBundleYsu(SeqMatcher sid_5p, SeqMatcher umi_5p, SeqMatcher sid_3p, SeqMatcher umi_3p,
                           std::optional<std::string> bait_5p, std::optional<std::string> bait_3p)
    : _sid_5p_matcher{std::move(sid_5p)},
      _umi_5p_matcher{std::move(umi_5p)},
      _sid_3p_matcher{std::move(sid_3p)},
      _umi_3p_matcher{std::move(umi_3p)},
      _bait_5p{std::move(bait_5p)},
      _bait_3p{std::move(bait_3p)} {}

const SeqMatcher& LutBundleYsu::Sid5pMatcher() const { return _sid_5p_matcher; }

const BarcodePool& LutBundleYsu::Sid5pPool() const { return _sid_5p_matcher.Pool(); }

const SeqMatcher& LutBundleYsu::Umi5pMatcher() const { return _umi_5p_matcher; }

const BarcodePool& LutBundleYsu::Umi5pPool() const { return _umi_5p_matcher.Pool(); }

const SeqMatcher& LutBundleYsu::Sid3pMatcher() const { return _sid_3p_matcher; }

const BarcodePool& LutBundleYsu::Sid3pPool() const { return _sid_3p_matcher.Pool(); }

const SeqMatcher& LutBundleYsu::Umi3pMatcher() const { return _umi_3p_matcher; }

const BarcodePool& LutBundleYsu::Umi3pPool() const { return _umi_3p_matcher.Pool(); }

const std::optional<std::string>& LutBundleYsu::Bait5p() const { return _bait_5p; }

const std::optional<std::string>& LutBundleYsu::Bait3p() const { return _bait_3p; }

template <>
LutBundleYsu CreateLutBundle(const AdapterDesignBundle& designs) {
  const auto& design = designs.design;
  if (design.bait_5p.has_value() && design.bait_5p->empty()) {
    throw error::Error("Adapter design '{}': bait_5p must not be empty", design.name);
  }
  if (design.bait_3p.has_value() && design.bait_3p->empty()) {
    throw error::Error("Adapter design '{}': bait_3p must not be empty", design.name);
  }
  return LutBundleYsu{
      *designs.GetMatcher5p(BarcodeType::kSid),
      *designs.GetMatcher5p(BarcodeType::kUmi),
      *designs.GetMatcher3p(BarcodeType::kSid),
      *designs.GetMatcher3p(BarcodeType::kUmi),
      design.bait_5p,
      design.bait_3p,
  };
}

}  // namespace xoos::demux
