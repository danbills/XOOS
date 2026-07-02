#pragma once

#include <optional>
#include <string>

#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

namespace fs = std::filesystem;

/**
 * A bundle of LUT to be used in trimming YSU data. These LUTs are meant
 * to be loaded from disk, but they might be created in memory for tests cases.
 */
class LutBundleYsu {
 public:
  LutBundleYsu(SeqMatcher sid_5p, SeqMatcher umi_5p, SeqMatcher sid_3p, SeqMatcher umi_3p,
               std::optional<std::string> bait_5p, std::optional<std::string> bait_3p);

  const SeqMatcher& Sid5pMatcher() const;
  const BarcodePool& Sid5pPool() const;

  const SeqMatcher& Umi5pMatcher() const;
  const BarcodePool& Umi5pPool() const;

  const SeqMatcher& Sid3pMatcher() const;
  const BarcodePool& Sid3pPool() const;

  const SeqMatcher& Umi3pMatcher() const;
  const BarcodePool& Umi3pPool() const;

  const std::optional<std::string>& Bait5p() const;
  const std::optional<std::string>& Bait3p() const;

 private:
  SeqMatcher _sid_5p_matcher;
  SeqMatcher _umi_5p_matcher;

  SeqMatcher _sid_3p_matcher;
  SeqMatcher _umi_3p_matcher;

  std::optional<std::string> _bait_5p;
  std::optional<std::string> _bait_3p;
};
}  // namespace xoos::demux
