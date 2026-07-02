#include "demux-and-trim-ysu.h"

#include <xoos/types/int.h>

#include <algorithm>

#include "io/read-record.h"

namespace xoos::demux {

static Trim5pYsu CreateTrim5p(const bool enable_partial, const LutBundleYsu& lut_bundle) {
  return Trim5pYsu{enable_partial, lut_bundle.Sid5pMatcher(), lut_bundle.Umi5pMatcher(), lut_bundle.Bait5p()};
}

static Trim3pYsu CreateTrim3p(const bool enable_partial, const LutBundleYsu& lut_bundle) {
  return Trim3pYsu{enable_partial, lut_bundle.Sid3pMatcher(), lut_bundle.Umi3pMatcher(), lut_bundle.Bait3p()};
}

DemuxAndTrimYsu::DemuxAndTrimYsu(const bool enable_partial, const LutBundleYsu& lut_bundle)
    : _enable_partial{enable_partial},
      _trim_5p{CreateTrim5p(enable_partial, lut_bundle)},
      _trim_3p{CreateTrim3p(enable_partial, lut_bundle)} {}

TrimInfoYsu DemuxAndTrimYsu::operator()(const FixedReadRecord& record) const {
  // speed optimization by converting the sequence to a 2-bit representation.
  const auto trim_5p = _trim_5p.Trim(record.TwoBitsSeq(), record.SeqLen(), record.Seq());
  const auto trim_3p = _trim_3p.Trim(record.TwoBitsSeq(), record.SeqLen(), record.Seq(), trim_5p.insert_start);
  return Demux(trim_5p, trim_3p, record.SeqLen());
}

TrimInfoYsu DemuxAndTrimYsu::Demux(const Trim5pInfoYsu& trim_5p, const Trim3pInfoYsu& trim_3p,
                                   const u32 seq_length) const {
  const auto sid_5p = MatchInfo::BarcodeId(trim_5p.sid_match);
  const auto sid_3p = MatchInfo::BarcodeId(trim_3p.sid_match);
  const auto sid_5p_edist = MatchInfo::EDist(trim_5p.sid_match);
  const auto sid_3p_edist = MatchInfo::EDist(trim_3p.sid_match);

  auto insert_start = trim_5p.insert_start;
  auto insert_end = trim_3p.insert_end;

  auto umi_5p = MatchInfo::BarcodeId(trim_5p.umi_match);
  auto umi_3p = MatchInfo::BarcodeId(trim_3p.umi_match);

  // Winning-side-only trimming for discordant reads: clear the losing side's trim and UMI
  if (sid_5p && sid_3p && *sid_5p != *sid_3p) {
    if (sid_5p_edist && sid_3p_edist && *sid_3p_edist < *sid_5p_edist) {
      // 3' wins — don't trim 5', clear 5' UMI
      insert_start = 0;
      umi_5p = std::nullopt;
    } else {
      // 5' wins (including tie) — don't trim 3', clear 3' UMI
      insert_end = seq_length;
      umi_3p = std::nullopt;
    }
  }

  return TrimInfoYsu{
      DetermineSampleId(trim_5p, trim_3p),
      sid_5p,
      sid_5p_edist,
      umi_5p,
      sid_3p,
      sid_3p_edist,
      umi_3p,
      LociRange{std::min(insert_start, insert_end), insert_end},
  };
}

std::optional<u32> DemuxAndTrimYsu::DetermineSampleId(const Trim5pInfoYsu& trim_5p,
                                                      const Trim3pInfoYsu& trim_3p) const {
  const auto has_umi =
      _enable_partial ? trim_5p.umi_match || trim_3p.umi_match : trim_5p.umi_match && trim_3p.umi_match;
  if (!has_umi) {
    return std::nullopt;
  }

  const auto sid_edist_5p = MatchInfo::EDist(trim_5p.sid_match);
  const auto sid_edist_3p = MatchInfo::EDist(trim_3p.sid_match);

  if (sid_edist_5p && !sid_edist_3p) {
    return MatchInfo::BarcodeId(trim_5p.sid_match);
  }

  if (!sid_edist_5p && sid_edist_3p) {
    return MatchInfo::BarcodeId(trim_3p.sid_match);
  }

  if (!sid_edist_5p && !sid_edist_3p) {
    return std::nullopt;
  }

  if (sid_edist_3p < sid_edist_5p) {
    return MatchInfo::BarcodeId(trim_3p.sid_match);
  }

  return MatchInfo::BarcodeId(trim_5p.sid_match);
}

}  // namespace xoos::demux
