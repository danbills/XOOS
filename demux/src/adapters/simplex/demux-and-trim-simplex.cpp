#include "demux-and-trim-simplex.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <string_view>
#include <tuple>

#include <gtl/phmap.hpp>
#include <magic_enum/magic_enum.hpp>

#include "xoos/log/logging.h"

namespace xoos::demux {

DemuxAndTrimSimplex::DemuxAndTrimSimplex(const DemuxAndTrimParam& param, const LutBundleSimplex& lut_bundle)
    : _param{param},
      _lut_bundle{lut_bundle},
      _nominal_seed_size{static_cast<s32>(lut_bundle.sid_2_matcher.Lut()->Pool().front().sequence.size())},
      _bw_extend_1_dfa{_lut_bundle.fixed_1_matcher.GetFrontSeq(), SearchDirection::kBackward},
      // TODO: change max_drop_difference hardcoding here
      _fw_extend_3_dfa{_lut_bundle.fixed_3_matcher.GetFrontSeq(), SearchDirection::kForward, 10},
      _bw_extend_5_dfa{_lut_bundle.fixed_5_matcher.GetFrontSeq(), SearchDirection::kBackward, 10},
      _fw_extend_7_dfa{_lut_bundle.fixed_7_matcher.GetFrontSeq()},
      _fw_search_1_bitap(_lut_bundle.fixed_1_matcher.GetFrontSeq(), SearchDirection::kForward),
      _bw_search_7_bitap(_lut_bundle.fixed_7_matcher.GetFrontSeq(), SearchDirection::kBackward),
      _sid_min_k(sid_scoring::SidScoring::CalcSidMinK(lut_bundle.sid_2_matcher)),
      _sid_max_k(
          static_cast<u8>(_lut_bundle.sid_2_matcher.GetFrontSeq().size() + _lut_bundle.sid_2_matcher.GetMaxDist())),
      _sid_kmer_prob(sid_scoring::SidScoring::CalcSIDProbRandomMatch(lut_bundle.sid_2_matcher, lut_bundle.sid_6_matcher,
                                                                     _sid_min_k, _sid_max_k)),
      _search_space_coefficients(
          sid_scoring::SidScoring::CalcSearchSpaceCoefficients(_sid_kmer_prob, _sid_min_k, _sid_max_k)) {
  // Validate that 5' and 3' SID pools have matching k-mer size ranges.
  sid_scoring::SidScoring::ValidateSidBounds(_lut_bundle.sid_6_matcher);
  const auto nominal_kmer_size = _lut_bundle.sid_6_matcher.GetFrontSeq().size();
  const auto edist = _lut_bundle.sid_6_matcher.GetMaxDist();

  const auto second_sid_min_len = nominal_kmer_size - edist;
  const auto second_sid_max_len = nominal_kmer_size + edist;

  // SID luts should have the same min and max_k
  if (_sid_min_k != second_sid_min_len) {
    throw error::Error("SIDs for simplex adapter need to be the min length. Saw {} and {} ", _sid_min_k,
                       second_sid_min_len);
  }
  if (_sid_max_k != second_sid_max_len) {
    throw error::Error("SIDs for simplex adapter need to be the max length. Saw {} and {} ", _sid_max_k,
                       second_sid_max_len);
  }

  const auto min_length_score_threshold =
      param.min_score +
      sid_scoring::SidScoring::CalcSearchSpacePenalty(_search_space_coefficients, static_cast<u32>(param.min_read_len));

  const auto min_match_length = static_cast<u32>((min_length_score_threshold + scoring::kMatch - 1) / scoring::kMatch);

  if (param.min_read_len < min_match_length) {
    throw error::Error(
        "The maximum seed sid search sequence size is {} the --min_read_length is {} which is pointless to process. "
        "set --min_read_length to a value >= {} for this adapter.",
        min_match_length, param.min_read_len, min_match_length);
  }

  Logging::Debug(
      "Running adapter design with min_score of {} and {} sids. This yields a score "
      "threshold of {} for sequences {} bp long and score of {} for sequences {} bp long.",
      param.min_score, _lut_bundle.sid_2_matcher.Pool().size(), min_length_score_threshold, param.min_read_len,
      param.min_score + sid_scoring::SidScoring::CalcSearchSpacePenalty(_search_space_coefficients, kMaxReadLength),
      kMaxReadLength);
}

/**
 * @brief Converts per-k unique k-mer counts into random match probabilities.
 *
 * @param unique_kmer_counts  Number of unique k-mers for each k in [min_k, max_k], indexed by (k - min_k).
 * @param min_k  Minimum k-mer size.
 * @param max_k  Maximum k-mer size.
 * @return Vector of probabilities indexed by (k − min_k).
 */
static std::vector<f64> KmerCountsToProb(const std::vector<f64>& unique_kmer_counts, const u8 min_k, const u8 max_k) {
  std::vector<f64> sid_kmer_prob(max_k - min_k + 1);
  for (u8 k = min_k; k <= max_k; ++k) {
    const auto unique_count = unique_kmer_counts[k - min_k];
    const auto kmer_count = std::ldexp(1.0, 2 * k);
    sid_kmer_prob[k - min_k] = unique_count / kmer_count;
    Logging::Debug(
        "For k-mer size {}, unique k-mers in SIDs: {}, total possible k-mers: {}, random match probability: {}", k,
        static_cast<u64>(unique_count), static_cast<u64>(kmer_count), sid_kmer_prob[k - min_k]);
  }
  return sid_kmer_prob;
}

std::vector<f64> sid_scoring::SidScoring::CalcSIDProbRandomMatch(const SeqMatcher& sid_matcher, const u8 min_k,
                                                                 const u8 max_k) {
  std::vector<f64> counts(max_k - min_k + 1);
  for (u8 k = min_k; k <= max_k; ++k) {
    counts[k - min_k] = static_cast<f64>(sid_matcher.Lut()->HashTables().at(k).size());
  }
  return KmerCountsToProb(counts, min_k, max_k);
}

std::vector<f64> sid_scoring::SidScoring::CalcSIDProbRandomMatch(const SeqMatcher& sid_matcher_a,
                                                                 const SeqMatcher& sid_matcher_b, const u8 min_k,
                                                                 const u8 max_k) {
  std::vector<f64> counts(max_k - min_k + 1);
  for (u8 k = min_k; k <= max_k; ++k) {
    counts[k - min_k] = static_cast<f64>(sid_matcher_a.Lut()->HashTables().at(k).size() +
                                         sid_matcher_b.Lut()->HashTables().at(k).size());
  }
  return KmerCountsToProb(counts, min_k, max_k);
}

sid_scoring::SidScoring::SearchSpaceCoefficients sid_scoring::SidScoring::CalcSearchSpaceCoefficients(
    const std::vector<f64>& sid_kmer_prob, const u8 min_k, const u8 max_k) {
  f64 linear = 0.0;
  f64 constant = 0.0;
  for (u8 k = min_k; k <= max_k; ++k) {
    const auto p = sid_kmer_prob[k - min_k];
    linear += p;
    constant += p * (1.0 - k);
  }
  Logging::Debug("Precomputed search-space coefficients: linear: {}, constant: {}", linear, constant);
  return {linear, constant};
}

s32 sid_scoring::SidScoring::CalcSearchSpacePenalty(const SearchSpaceCoefficients& coefficients, const u32 seq_len) {
  const auto expected_random_hits = coefficients.linear * seq_len + coefficients.constant;
  return static_cast<s32>(std::ceil(std::log2(expected_random_hits > 1.0 ? expected_random_hits : 1.0)));
}

void sid_scoring::SidScoring::ValidateSidBounds(const SeqMatcher& sid_matcher) {
  const auto sid_len = sid_matcher.GetFrontSeq().size();
  const auto sid_edist = sid_matcher.GetMaxDist();
  if (sid_len < sid_edist) {
    throw error::Error("SID length {} is shorter than edit distance {}", sid_len, sid_edist);
  }
}

u8 sid_scoring::SidScoring::CalcSidMinK(const SeqMatcher& sid_matcher) {
  ValidateSidBounds(sid_matcher);
  return static_cast<u8>(sid_matcher.GetFrontSeq().size() - sid_matcher.GetMaxDist());
}

const TrimInfoSimplex& DemuxAndTrimSimplex::operator()(FixedReadRecord& record) const {
  const auto length = record.SeqLen();

  // The surviving insert of our search [spos, epos) must be at least `min_insert_len` bases long (matching the
  // downstream filter in demux.cpp). Stated once here so both ends share it and it cannot drift or be inverted.
  const auto min_insert_len = static_cast<s32>(_param.min_trimmed_read_len);

  // checks if our insert is long enough if not we can reject additional searches
  const auto leaves_enough_insert = [min_insert_len](const s32 spos, const s32 epos) {
    return epos - spos >= min_insert_len;
  };

  // deepest possible valid 5' trim
  const auto search_limit_5p = static_cast<s32>(length) - min_insert_len;
  // shallowest possible valid 3' trim
  const auto search_limit_3p = min_insert_len;

  const auto min_score =
      _param.min_score + sid_scoring::SidScoring::CalcSearchSpacePenalty(_search_space_coefficients, length);

  using SID = u32;
  using enum DFAClassifier::DFAResult::State;
  // barcode_id -> (MatchInfo, DFAResults, score, trim_pos)
  gtl::flat_hash_map<SID, std::tuple<MatchInfo, PairedDFAResults, s32, s32>> candidates_5p;
  candidates_5p.reserve(_lut_bundle.sid_2_matcher.Pool().size());

  auto best_5p_score = std::numeric_limits<s32>::min();
  constexpr SID kNoSID = std::numeric_limits<SID>::max();
  SID best_5p_sid = kNoSID;
  // Step 1: Search from the 5' side for SID_2 to seed the match until
  {
    // TODO: Can perform baiting to find candidate search position to prevent unneeded work (Speed improvement)
    auto search_pos = 0;
    std::vector<MatchInfo> sid_matches;

    while (search_pos < search_limit_5p) {
      // greedy search for next barcode to see match
      _lut_bundle.sid_2_matcher.FindNextBarcode(ReadEnd::k5p, search_pos, record.TwoBitsSeq(),
                                                static_cast<size_t>(search_limit_5p), sid_matches);
      if (sid_matches.empty()) {
        // finished processing as we don't can't find a match
        break;
      }
      search_pos = static_cast<s32>(sid_matches[0].SPos()) + 1;
      for (auto sid_match : sid_matches) {
        // Extend with backward DFA (fixed_1 region, left of SID) and forward DFA (fixed_3 region, right of SID)
        auto dfa_results = ExtendCandidate(record, sid_match, _bw_extend_1_dfa, _fw_extend_3_dfa);
        auto score = sid_match.LogLikelihoodScore(_nominal_seed_size) + dfa_results.bw_dfa_results.best.score;
        auto trim_pos_5p = static_cast<s32>(sid_match.EPos());

        // calc with last_match towards insert
        if (dfa_results.fw_dfa_results.last_match.has_value()) {
          const auto& last = dfa_results.fw_dfa_results.last_match.value();
          // add in the extra bases that we matched
          trim_pos_5p += last.pos;
          switch (last.state) {
            case kFoundEnd:
              // we found the end of the dfa match
              score += last.score;
              break;
            case kUnresolved: {
              const auto remaining_bases = dfa_results.fw_dfa_results.remaining_last_match_bases;
              // add in remaining bases as mismatches
              score += last.score + remaining_bases * scoring::kSubstitution;
              trim_pos_5p += remaining_bases;
              break;
            }
            case kFoundPruned:
            case kInputEnd: {
              // somehow we made it off the end (likely truncated) and the remaining bases must be deletions
              const auto remaining_bases = dfa_results.fw_dfa_results.remaining_last_match_bases;
              score += last.score + remaining_bases * scoring::kDeletion;
              break;
            }
            default:
              throw error::Error("Found impossible DFA state end {} ", magic_enum::enum_name(last.state));
          }
        } else {
          const auto remaining_bases = dfa_results.fw_dfa_results.remaining_last_match_bases;
          // add in remaining bases as mismatches
          score += remaining_bases * scoring::kSubstitution;
          trim_pos_5p += remaining_bases;
        }

        if (!leaves_enough_insert(trim_pos_5p, static_cast<s32>(length))) {
          // 3' boundary is still the read end here; trimming this far would leave too little insert.
          continue;
        }

        // insert new entry or update existing if this score is higher
        const auto [it, inserted] =
            candidates_5p.try_emplace(sid_match.BarcodeId(), sid_match, dfa_results, score, trim_pos_5p);
        if (!inserted && score > std::get<2>(it->second)) {
          it->second = {sid_match, dfa_results, score, trim_pos_5p};
        }

        if (score > best_5p_score) {
          best_5p_score = score;
          best_5p_sid = sid_match.BarcodeId();
          // TODO: Can threshold for speed here if we have a clear match without random chance (Speed improvement)
          // if (score >= kScoreThreshold) {
          //   score_threshold_reached = true;
          // }
        }
      }
    }
  }
  // test singleton 3p candidate
  auto best_3p_score = std::numeric_limits<s32>::min();
  MatchInfo best_3p_match;
  auto best_3p_trim_pos = 0;

  // best pair candidate
  auto& best_pair_candidate = record.trim_info_simplex;
  // reset loci
  best_pair_candidate.insert.epos = length;

  auto best_pair_score = std::numeric_limits<s32>::min();

  // Step 2: Search for 3' SID_6 from the 3' side
  {
    auto search_pos = static_cast<s32>(length);
    std::vector<MatchInfo> sid_matches;
    while (search_pos > search_limit_3p) {
      _lut_bundle.sid_6_matcher.FindNextBarcode(ReadEnd::k3p, search_pos, record.TwoBitsSeq(), length, sid_matches);
      if (sid_matches.empty()) {
        // finished processing as we don't can't find a match
        break;
      }
      search_pos = static_cast<s32>(sid_matches[0].EPos()) - 1;
      for (const auto sid_match : sid_matches) {
        // Extend with backward DFA (fixed_1 region, left of SID) and forward DFA (fixed_3 region, right of SID)
        const auto dfa_results = ExtendCandidate(record, sid_match, _bw_extend_5_dfa, _fw_extend_7_dfa);
        auto score = sid_match.LogLikelihoodScore(_nominal_seed_size) + dfa_results.fw_dfa_results.best.score;
        auto trim_pos_3p = static_cast<s32>(sid_match.SPos());

        // calc with last_match towards insert
        if (dfa_results.bw_dfa_results.last_match.has_value()) {
          const auto& last = dfa_results.bw_dfa_results.last_match.value();
          // add in the extra bases that we matched
          trim_pos_3p -= last.pos;
          switch (last.state) {
            case kFoundEnd:
              // we found the end of the dfa match
              score += last.score;
              break;
            case kInputEnd: {
              // somehow we made it off the end (likely truncated) and the remaining bases must be deletions
              const auto remaining_bases = dfa_results.bw_dfa_results.remaining_last_match_bases;
              score += last.score + remaining_bases * scoring::kDeletion;
              break;
            }
            case kFoundPruned:
            case kUnresolved: {
              // we didn't find any meaningful end
              const auto remaining_bases = dfa_results.bw_dfa_results.remaining_last_match_bases;
              // add in remaining bases as mismatches
              score += last.score + remaining_bases * scoring::kSubstitution;
              trim_pos_3p -= remaining_bases;
              break;
            }
            default:
              throw error::Error("Trim location ambiguous. Found state end {} ", magic_enum::enum_name(last.state));
          }
        } else {
          const auto remaining_bases = dfa_results.bw_dfa_results.remaining_last_match_bases;
          // add in remaining bases as mismatches
          score += remaining_bases * scoring::kSubstitution;
          trim_pos_3p -= remaining_bases;
        }

        if (!leaves_enough_insert(0, trim_pos_3p)) {
          // 5' boundary is still 0 here; trimming this far would leave too little insert.
          continue;
        }

        const auto match_5p = candidates_5p.find(sid_match.BarcodeId());

        // if we find a paired SID assign it
        if (match_5p != candidates_5p.end()) {
          const auto& [match_info_5p, dfa_results_5p, score_5p, trim_pos_5p] = match_5p->second;
          const auto pair_score = score + score_5p;
          // pair candidate only added if the join score is higher than threshold
          if (pair_score > best_pair_score && trim_pos_5p < trim_pos_3p && pair_score >= min_score) {
            best_pair_candidate.sid = sid_match.BarcodeId();
            best_pair_candidate.sid_5p = match_info_5p.BarcodeId();
            best_pair_candidate.sid_5p_edist = match_info_5p.EDist();
            best_pair_candidate.score_5p = score_5p;
            best_pair_candidate.sid_3p = sid_match.BarcodeId();
            best_pair_candidate.sid_3p_edist = sid_match.EDist();
            best_pair_candidate.score_3p = score;
            best_pair_candidate.insert.spos = static_cast<u32>(trim_pos_5p);
            best_pair_candidate.insert.epos = static_cast<u32>(trim_pos_3p);
            best_pair_score = pair_score;
          }
        } else if (best_3p_score < score) {
          best_3p_score = score;
          best_3p_match = sid_match;
          best_3p_trim_pos = trim_pos_3p;
          // TODO: Can threshold for speed here if we have a clear match without random chance (Speed improvement)
        }
      }
    }
  }

  // Stage 3 consider the best candidates to return
  if (best_pair_score >= std::max(best_5p_score, best_3p_score)) {
    // full read that no partial beats
    return best_pair_candidate;
  }

  // clear so we can reconstruct a partial read
  best_pair_candidate.Clear();

  // if the paired event is worse than a singleton score
  // Check if partial read is possible
  if (best_5p_score < min_score && best_3p_score < min_score) {
    // return empty candidate
    return best_pair_candidate;
  }
  if (best_5p_score >= min_score) {
    const auto& [match_info_5p, dfa_results_5p, score_5p, trim_pos_5p] = candidates_5p.at(best_5p_sid);
    best_pair_candidate.sid_5p = match_info_5p.BarcodeId();
    best_pair_candidate.sid_5p_edist = match_info_5p.EDist();
    best_pair_candidate.score_5p = best_5p_score;
  }
  if (best_3p_score >= min_score) {
    best_pair_candidate.sid_3p = best_3p_match.BarcodeId();
    best_pair_candidate.sid_3p_edist = best_3p_match.EDist();
    best_pair_candidate.score_3p = best_3p_score;
  }
  // if we allow partial reads then we can assign an SID and its trimming position
  if (_param.read_length_mode != ReadLengthMode::kFullOnly) {
    if (best_3p_score > best_5p_score) {
      // 3' wins: trim 3' side, leave 5' untrimmed
      best_pair_candidate.sid = best_3p_match.BarcodeId();
      best_pair_candidate.insert.spos = 0;
      best_pair_candidate.insert.epos = static_cast<u32>(best_3p_trim_pos);
    } else {
      // 5' wins (including tie): trim 5' side, leave 3' untrimmed
      const auto& [match_info_5p, dfa_results_5p, score_5p, trim_pos_5p] = candidates_5p.at(best_5p_sid);
      best_pair_candidate.sid = best_5p_sid;
      best_pair_candidate.insert.spos = static_cast<u32>(trim_pos_5p);
      best_pair_candidate.insert.epos = length;
    }
  }
  return best_pair_candidate;
}

DemuxAndTrimSimplex::PairedDFAResults DemuxAndTrimSimplex::ExtendCandidate(const FixedReadRecord& record,
                                                                           const MatchInfo& match_info,
                                                                           const DFAClassifier& bw_dfa,
                                                                           const DFAClassifier& fw_dfa) {
  const auto search_domain_bw = std::string_view(record.Seq(), match_info.SPos());
  const auto search_domain_fw = std::string_view(record.Seq() + match_info.EPos(), record.SeqLen() - match_info.EPos());
  return {bw_dfa.Classify(search_domain_bw), fw_dfa.Classify(search_domain_fw)};
}

}  // namespace xoos::demux
