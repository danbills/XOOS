#pragma once

#include <optional>
#include <string>

#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

namespace fs = std::filesystem;

/**
 * A bundle of LUT to be used in trimming YS data. These LUTs are meant
 * to be loaded from disk, but they might be created in memory for tests cases.
 */
class LutBundleYs {
 public:
  LutBundleYs(SeqMatcher runway_5p, SeqMatcher sid_5p, SeqMatcher sid_spacer_5p, SeqMatcher sid_3p,
              std::optional<std::string> bait_3p);

  const SeqMatcher& Runway5pMatcher() const;

  const SeqMatcher& Sid5pMatcher() const;
  const BarcodePool& Sid5pPool() const;

  const SeqMatcher& SidSpacer5pMatcher() const;
  const BarcodePool& SidSpacer5pPool() const;

  const SeqMatcher& Sid3pMatcher() const;
  const BarcodePool& Sid3pPool() const;

  const std::optional<std::string>& Bait3p() const;

 private:
  SeqMatcher _runway_5p_matcher;
  SeqMatcher _sid_5p_matcher;
  SeqMatcher _sid_spacer_5p_matcher;

  SeqMatcher _sid_3p_matcher;

  std::optional<std::string> _bait_3p;
};
}  // namespace xoos::demux
