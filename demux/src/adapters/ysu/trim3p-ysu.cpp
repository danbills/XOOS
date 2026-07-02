#include "trim3p-ysu.h"

#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>

#include <algorithm>
#include <string_view>
#include <utility>

namespace xoos::demux {

Trim3pInfoYsu::Trim3pInfoYsu(const u32 insert_end) : insert_end(insert_end) {}

Trim3pYsu::Trim3pYsu(const bool enable_partial, SeqMatcher sid_3p_matcher, SeqMatcher umi_3p_matcher,
                     std::optional<std::string> bait_3p)
    : _enable_partial(enable_partial),
      _sid_3p_matcher(std::move(sid_3p_matcher)),
      _umi_3p_matcher(std::move(umi_3p_matcher)),
      _bait_3p{std::move(bait_3p)} {}

Trim3pInfoYsu Trim3pYsu::Trim(const u8* const seq2, const u32 seq_length, const char* const seq,
                              const u32 insert_start) const {
  Trim3pInfoYsu trim_info(seq_length);
  std::optional<MatchInfo> sid_match;
  u32 current_pos = seq_length;

  // Search for the bait sequence adjacent to the end of 3'sid and update
  // current position. This would better handle situations where there are more
  // extra bases after the SID than expected.
  if (_bait_3p.has_value()) {
    const u32 search_length = std::min(64U, seq_length);
    const std::string_view tail(seq + seq_length - search_length, search_length);
    const auto pos = tail.rfind(_bait_3p.value());
    if (pos != std::string_view::npos) {
      current_pos = seq_length - search_length + static_cast<u32>(pos);
    }
  }

  {
    auto sid_match0 = _sid_3p_matcher.FindBarcode(ReadEnd::k3p, current_pos, seq2, seq_length);
    if (AbortTrim(sid_match0, insert_start)) {
      return trim_info;
    }

    if (!sid_match0.IsUnknown()) {
      current_pos = sid_match0.SPos();
      sid_match = sid_match0;
    }
  }

  {
    auto umi_match = _umi_3p_matcher.FindBarcode(ReadEnd::k3p, current_pos, seq2, seq_length);
    if (AbortTrim(umi_match, insert_start)) {
      return trim_info;
    }

    // Only report SID for trimming when UMI is also found on this end.
    // A SID-only match is unreliable for trimming because the 12bp SID
    // can match spuriously within the genomic insert.
    if (!umi_match.IsUnknown()) {
      current_pos = umi_match.SPos();
      trim_info.sid_match = sid_match;
      trim_info.umi_match = umi_match;
      trim_info.insert_end = current_pos;
    }
  }

  return trim_info;
}

bool Trim3pYsu::AbortTrim(const MatchInfo& match_info, const u32 insert_start) const {
  return (!match_info.IsUnknown() && match_info.SPos() <= insert_start) ||
         (match_info.IsUnknown() && !_enable_partial) || match_info.IsAmbiguous();
}

}  // namespace xoos::demux
