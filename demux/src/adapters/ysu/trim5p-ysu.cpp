#include "trim5p-ysu.h"

#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>

#include <algorithm>
#include <string_view>
#include <utility>

namespace xoos::demux {

Trim5pYsu::Trim5pYsu(const bool enable_partial, SeqMatcher sid_5p_matcher, SeqMatcher umi_5p_matcher,
                     std::optional<std::string> bait_5p)
    : _enable_partial{enable_partial},
      _sid_5p_matcher(std::move(sid_5p_matcher)),
      _umi_5p_matcher(std::move(umi_5p_matcher)),
      _bait_5p{std::move(bait_5p)} {}

Trim5pInfoYsu Trim5pYsu::Trim(const u8* const seq2, const u32 length, const char* const seq) const {
  u32 current_pos{0};

  Trim5pInfoYsu trim_info;

  // Search for the bait sequence adjacent to the start of 5'sid and update
  // current position. This would better handle situations where there are more
  // extra bases before the SID than expected.
  if (_bait_5p.has_value()) {
    const u32 search_length = std::min(64U, length);
    const std::string_view head(seq, search_length);
    const auto pos = head.find(_bait_5p.value());
    if (pos != std::string_view::npos) {
      current_pos = static_cast<u32>(pos) + static_cast<u32>(_bait_5p->size());
    }
  }

  std::optional<MatchInfo> sid_match;
  {
    auto sid_match0 = _sid_5p_matcher.FindBarcode(ReadEnd::k5p, current_pos, seq2, length);
    if (AbortTrim(sid_match0)) {
      return trim_info;
    }

    if (!sid_match0.IsUnknown()) {
      current_pos = sid_match0.EPos();
      sid_match = sid_match0;
    }
  }

  {
    auto umi_match = _umi_5p_matcher.FindBarcode(ReadEnd::k5p, current_pos, seq2, length);
    if (AbortTrim(umi_match)) {
      return trim_info;
    }

    // Only report SID for trimming when UMI is also found on this end.
    // A SID-only match is unreliable for trimming because the 12bp SID
    // can match spuriously within the genomic insert.
    if (!umi_match.IsUnknown()) {
      current_pos = umi_match.EPos();
      trim_info.sid_match = sid_match;
      trim_info.umi_match = umi_match;
      trim_info.insert_start = current_pos;
    }
  }

  return trim_info;
}

bool Trim5pYsu::AbortTrim(const MatchInfo& match_info) const {
  return (match_info.IsUnknown() && !_enable_partial) || match_info.IsAmbiguous();
}
}  // namespace xoos::demux
