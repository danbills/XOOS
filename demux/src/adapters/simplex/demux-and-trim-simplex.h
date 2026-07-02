#pragma once

#include "core/demux-and-trim-pipeline.h"
#include "io/read-record.h"
#include "lut-bundle-simplex.h"
#include "sequence/matcher/bitap.h"
#include "sequence/matcher/dfa-classifier.h"
#include "trim-info-simplex.h"

namespace xoos::demux {

namespace sid_scoring {

/**
 * @class SidScoring
 * @brief Shared SID scoring utilities for search-space penalty calculation.
 *
 * Provides Karlin-Altschul-style search-space correction for LUT-based SID matching.
 * Used by both simplex and simplex-10x adapters.
 */
class SidScoring {
 public:
  /// Precomputed linear regression coefficients for O(1) search-space penalty calculation.
  struct SearchSpaceCoefficients {
    f64 linear;
    f64 constant;
  };

  /**
   * @brief Computes per-k probability of a random k-mer matching any SID in the LUT.
   *
   * For each k in [min_k, max_k], calculates unique_kmers / 4^k using the SID LUT hash table.
   *
   * @param sid_matcher  The SID matcher whose LUT hash tables are queried.
   * @param min_k        Minimum k-mer size to evaluate.
   * @param max_k        Maximum k-mer size to evaluate.
   * @return Vector of probabilities indexed by (k − min_k).
   */
  static std::vector<f64> CalcSIDProbRandomMatch(const SeqMatcher& sid_matcher, u8 min_k, u8 max_k);

  /**
   * @brief Computes per-k probability of a random k-mer matching any SID in two LUTs.
   *
   * For each k in [min_k, max_k], calculates (count_a + count_b) / 4^k using both SID LUT hash
   * tables. The sum is a conservative upper bound (not a union) that may slightly overestimate
   * when k-mers overlap between the two pools.
   * Used by adapters with separate 5' and 3' SID pools (e.g., simplex/YS-NEW).
   *
   * @param sid_matcher_a  First SID matcher (e.g., sid_2).
   * @param sid_matcher_b  Second SID matcher (e.g., sid_6).
   * @param min_k          Minimum k-mer size to evaluate.
   * @param max_k          Maximum k-mer size to evaluate.
   * @return Vector of probabilities indexed by (k − min_k).
   */
  static std::vector<f64> CalcSIDProbRandomMatch(const SeqMatcher& sid_matcher_a, const SeqMatcher& sid_matcher_b,
                                                 u8 min_k, u8 max_k);

  /**
   * @brief Precomputes the linear and constant coefficients from per-k random match probabilities.
   *
   * These coefficients allow CalcSearchSpacePenalty to compute the expected random hit count in O(1),
   * analogous to precomputing the effective search-space parameters in Karlin-Altschul statistics.
   *
   * @param sid_kmer_prob  Per-k random match probabilities from CalcSIDProbRandomMatch.
   * @param min_k          Minimum k-mer size.
   * @param max_k          Maximum k-mer size.
   * @return SearchSpaceCoefficients with the precomputed linear and constant terms.
   */
  static SearchSpaceCoefficients CalcSearchSpaceCoefficients(const std::vector<f64>& sid_kmer_prob, u8 min_k, u8 max_k);

  /**
   * @brief Computes a length-dependent bit-score penalty to maintain consistent FDR per read.
   *
   * Conceptually equivalent to a Karlin-Altschul E-value correction: returns
   * ceil(log2(expected_random_hits)) for the given sequence length, i.e. the additional
   * base-2 bit-score required to offset the larger effective search space as read length grows.
   *
   * @param coefficients  Precomputed search-space coefficients from CalcSearchSpaceCoefficients.
   * @param seq_len       Length of the read sequence in bases.
   * @return The bit-score penalty (rounded up), always >= 0.
   */
  static s32 CalcSearchSpacePenalty(const SearchSpaceCoefficients& coefficients, u32 seq_len);

  /**
   * @brief Validates that a SID matcher's sequence length >= its max edit distance.
   *
   * Must be called before computing sid_min_k = len - edist to prevent unsigned underflow.
   *
   * @param sid_matcher  The SID matcher to validate.
   * @throws xoos::Error if the sequence length is shorter than the edit distance.
   */
  static void ValidateSidBounds(const SeqMatcher& sid_matcher);

  /**
   * @brief Computes sid_min_k = len - edist after validating bounds.
   *
   * @param sid_matcher  The SID matcher to compute min_k for.
   * @return The minimum k-mer size (sequence length minus max edit distance).
   */
  static u8 CalcSidMinK(const SeqMatcher& sid_matcher);
};

}  // namespace sid_scoring

/**
 * @class DemuxAndTrimSimplex
 * @brief Demultiplexes and trims simplex adapter reads.
 *
 * Identifies barcode SIDs at both the 5' and 3' ends of a read using LUT seed matching
 * followed by DFA-based extension. Computes alignment-like scores for each candidate,
 * resolves concordant and discordant pairs, and determines the insert (trimmed) region.
 */
class DemuxAndTrimSimplex {
 public:
  /**
   * @brief Constructs the demux-and-trim engine from runtime parameters and prebuilt LUT matchers.
   *
   * Initialises DFA classifiers, Bitap searchers, and precomputes score-modifier coefficients.
   * Validates that the SID lengths from both LUT bundles are consistent and that the minimum
   * read length is compatible with the scoring threshold.
   *
   * @param param  Runtime demux-and-trim parameters (score thresholds, read-length limits, etc.).
   * @param lut_bundle  Prebuilt LUT matchers and fixed-region sequences for the simplex adapter.
   * @throws error::Error If SID lengths are inconsistent or the score threshold is unreachable.
   */
  DemuxAndTrimSimplex(const DemuxAndTrimParam& param, const LutBundleSimplex& lut_bundle);

  /**
   * @brief Demultiplexes and trims a single read.
   *
   * Searches for 5' (SID_2) and 3' (SID_6) barcodes, extends each match with DFA classifiers,
   * scores them via log-likelihood alignment, and selects the best concordant pair or singleton.
   * The resulting trim coordinates and barcode IDs are written into the record's trim_info_simplex.
   *
   * @param record  The read record to process; trim_info_simplex is updated in place.
   * @return A const reference to the updated TrimInfoSimplex stored in the record.
   */
  const TrimInfoSimplex& operator()(FixedReadRecord& record) const;

 private:
  const DemuxAndTrimParam& _param;
  const LutBundleSimplex& _lut_bundle;

  /// @brief Nominal barcode seed length (bases), used for log-likelihood scoring.
  const s32 _nominal_seed_size;

  // --- DFA classifiers for extending LUT seed matches into flanking fixed regions ---
  /// @brief Backward DFA for extending into fixed region 1 (left of SID_2, 5' side).
  const DFAClassifier _bw_extend_1_dfa;
  /// @brief Forward DFA for extending into fixed region 3 (right of SID_2, towards insert).
  const DFAClassifier _fw_extend_3_dfa;
  /// @brief Backward DFA for extending into fixed region 5 (left of SID_6, towards insert).
  const DFAClassifier _bw_extend_5_dfa;
  /// @brief Forward DFA for extending into fixed region 7 (right of SID_6, 3' side).
  const DFAClassifier _fw_extend_7_dfa;

  // --- Bitap searchers for baiting faster complete searches ---
  /// @brief Forward Bitap searcher for fixed region 1.
  const Bitap<4> _fw_search_1_bitap;
  /// @brief Backward Bitap searcher for fixed region 7.
  const Bitap<4> _bw_search_7_bitap;

  /**
   * @brief Combined backward and forward DFA classification results for a single candidate.
   */
  struct PairedDFAResults {
    DFAClassifier::DFAResults bw_dfa_results;
    DFAClassifier::DFAResults fw_dfa_results;
  };

  /// @brief Minimum SID k-mer length (nominal length − max edit distance).
  const u8 _sid_min_k;
  /// @brief Maximum SID k-mer length (nominal length + max edit distance).
  const u8 _sid_max_k;
  /// @brief Per-k probability of a random k-mer matching any SID, indexed by (k − _sid_min_k).
  const std::vector<f64> _sid_kmer_prob;

  /**
   * @brief Precomputed coefficients for the O(1) effective search-space calculation
   *        in CalcSearchSpacePenalty.
   *
   * Analogous to the effective search-space term (m·n) in Karlin-Altschul statistics
   * (BLAST), the expected number of random k-mer hits for a sequence of length L is:
   *   expected_random_hits = linear * L + constant
   *
   * where:
   *   linear   = Σ_k prob[k]
   *   constant = Σ_k prob[k] * (1 - k)   (edge-effect correction, cf. effective length)
   *
   * Derived from expanding Σ_k prob[k] * (L - k + 1).
   */
  const sid_scoring::SidScoring::SearchSpaceCoefficients _search_space_coefficients;

  /**
   * @brief Extends a SID seed match in both directions using DFA classifiers.
   *
   * Runs the backward DFA over the region left of the seed and the forward DFA over the
   * region right of the seed to score how well the flanking fixed adapter regions match.
   *
   * @param record      The read being processed.
   * @param match_info  The SID seed match to extend.
   * @param bw_dfa      DFA classifier for the backward (left-of-seed) direction.
   * @param fw_dfa      DFA classifier for the forward (right-of-seed) direction.
   * @return Paired DFA classification results for both directions.
   */
  static PairedDFAResults ExtendCandidate(const FixedReadRecord& record, const MatchInfo& match_info,
                                          const DFAClassifier& bw_dfa, const DFAClassifier& fw_dfa);
};
}  // namespace xoos::demux
