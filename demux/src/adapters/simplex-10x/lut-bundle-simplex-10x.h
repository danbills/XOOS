#pragma once

#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

/**
 * LUT bundle for simplex-10x adapter (5'-only).
 *
 * Read structure: [fixed_1 (runway)][sid_2 (SID)][fixed_3 (stem)][insert...]
 * No 3' adapter components. The stem (fixed_3 / p_read_1t) is used for
 * scoring only — the trim position is set after the SID, not after the stem.
 */
class LutBundleSimplex10x {
 public:
  LutBundleSimplex10x(SeqMatcher fixed_1, SeqMatcher sid_2, SeqMatcher fixed_3);
  LutBundleSimplex10x(LutBundleSimplex10x&&) noexcept = default;
  LutBundleSimplex10x(const LutBundleSimplex10x&) = delete;
  LutBundleSimplex10x& operator=(LutBundleSimplex10x&&) = delete;
  LutBundleSimplex10x& operator=(const LutBundleSimplex10x&) = delete;

 public:
  // 1: anchor / runway (e.g. CAACAA)
  const SeqMatcher fixed_1_matcher;
  // 2: SID barcode (12 bp)
  const SeqMatcher sid_2_matcher;
  // 3: stem / p_read_1t (e.g. CTACACGACACTCTTCCGATCT, 22 bp) — scoring only, not trimmed.
  const SeqMatcher fixed_3_matcher;
};
}  // namespace xoos::demux
