#include "trim5p-ys.h"

#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>

#include <algorithm>
#include <string_view>
#include <utility>

#include "sequence/matcher/bitap.h"

namespace xoos::demux {

Trim5pYs::Trim5pYs(bool enable_partial, SeqMatcher runway_5p_matcher, SeqMatcher sid_5p_matcher,
                   SeqMatcher sid_spacer_5p_matcher)
    : _enable_partial{enable_partial},
      _runway_5p_matcher{std::move(runway_5p_matcher)},
      _sid_5p_matcher(std::move(sid_5p_matcher)),
      _sid_spacer_5p_matcher{std::move(sid_spacer_5p_matcher)},
      _flank_5p_sequence(_runway_5p_matcher.Pool().front().sequence + _sid_spacer_5p_matcher.Pool().front().sequence),
      _flank_5p_bitap(_flank_5p_sequence, SearchDirection::kForward) {}

/**
 * This is the main function that should be used to trim 5' adapters. It uses a 2-bit encoding (seq2) to allow
 * for fast LUT operations; it is using the original 8-bit encoding too for an exact string match.
 *
 * TODO: Why are we allowing for so much hardcoded strings that aren't using the bundle
 *       Alternatively, we could hardcode (if the optimizer can use this faster) but this should be via code generation
 *
 * @param seq2 2-bit encoded sequence for fast LUT operations
 * @param length length of the sequence
 * @param seq 8-bit encoded original sequence for exact string matching
 * @return Trim5pInfoYs structure containing trim information including SID match and insert start position
 */
Trim5pInfoYs Trim5pYs::Trim(const u8* seq2, size_t length, const char* seq) const {
  u32 current_pos{0};

  Trim5pInfoYs trim_info;

  if (length > 0) {
    const s32 end = std::min(Bitap<4>::kQueryWindowSize - 1, static_cast<s32>(length - 1));
    const s32 match_pos = _flank_5p_bitap.Find(seq, 0, end);
    if (match_pos != -1) {
      // Forward Bitap returns match end position; +1 gives position after flank.
      current_pos = static_cast<u32>(match_pos + 1);
    }
  }

  // Find the sid
  auto sid_match0 = _sid_5p_matcher.FindBarcode(ReadEnd::k5p, current_pos, seq2, length);
  if (AbortTrim(sid_match0)) {
    return trim_info;
  }

  if (!sid_match0.IsUnknown()) {
    current_pos = sid_match0.EPos();
    const std::optional sid_match = sid_match0;
    trim_info.sid_match = sid_match;
    trim_info.insert_start = current_pos;
  }
  return trim_info;
}

bool Trim5pYs::AbortTrim(const MatchInfo& match_info) const {
  return (match_info.IsUnknown() && !_enable_partial) || match_info.IsAmbiguous();
}
}  // namespace xoos::demux
