#pragma once

#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

/**
 * A generic simplex adapter originally designed for YS adapter but if similar in structure could be used for others.
 *
 * Can be initialized with any adapter with the structure with element in order:
 * [optional_0][fixed_1][sid_2][fixed_3][insert4][fixed_5][sid_6][fixed_7][optional_8]
 *
 * For example (lower case represents optional adapter sequence that may exist):
 * ttcagacgtgtactcttccgatctCAACAAGAGTCTTTT[sid_5p]GCGT[insert]CGC[sid_3p]TTTGTCGTGTAGGgaagagcacatgcattcgatgcggt
 *
 * In this example mapping:
 * optional_0 = ttcagacgtgtactcttccgatct
 * fixed_1 = CAACAAGAGTCTTTT
 * sid_2 = [sid_5p]
 * fixed_3 = GCGT
 * insert_4 = [insert]
 * fixed_5 = CGC
 * sid_6 = [sid_3p]
 * fixed_7 = TTTGTCGTGTAGG
 * optional_8 = gaagagcacatgcattcgatgcggt
 *
 * Currently we don't have objects representing insert or optional, though optional may be useful in the future.
 * The elements will be obtained upstream. Method doesn't matter but could be via adapter design bundle parsing.
 */
class LutBundleSimplex {
 public:
  LutBundleSimplex(SeqMatcher fixed_1, SeqMatcher sid_2, SeqMatcher fixed_3, SeqMatcher fixed_5, SeqMatcher sid_6,
                   SeqMatcher fixed_7);
  LutBundleSimplex(const LutBundleSimplex&) = default;
  LutBundleSimplex(LutBundleSimplex&&) noexcept = default;
  LutBundleSimplex& operator=(const LutBundleSimplex&) = delete;
  LutBundleSimplex& operator=(LutBundleSimplex&&) = delete;
  ~LutBundleSimplex() = default;

 public:
  // TODO: We don't actually use the fixed matchers, we only care about the fixed string derived from them via
  //       variable_name.Pool().front().sequence which is something we should change upstream
  // Variable are listed in order of the adapter components here:
  // 0 There may be unintended upstream adapter retained e.g. ttcagacgtgtactcttccgatct (E oligo), no object included yet
  // 1 fixed adapter should be here
  SeqMatcher fixed_1_matcher;
  // 2 variable adapter should be here
  SeqMatcher sid_2_matcher;
  // 3 fixed adapter should be here
  SeqMatcher fixed_3_matcher;
  // 4 Insert would be here, no object needed
  // 5 fixed adapter should be here
  SeqMatcher fixed_5_matcher;
  // 6 variable adapter should be here
  SeqMatcher sid_6_matcher;
  // 7 fixed adapter should be here
  SeqMatcher fixed_7_matcher;
  // 8 There may be unintended down adapter retained e.g. gaagagcacatgcattcgatgcggt (blocker), no object included yet
};
}  // namespace xoos::demux
