#include "intervals/interval-tools.h"

#include <algorithm>
#include <cmath>
#include <regex>

#include <fmt/core.h>
#include <utility/utility-functions.h>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/util/math.h>

namespace xoos::cnc {

/**
 * @brief merge intervals that are immediately adjacent to one another (but not overlapping). So (1,5) and (5,10) would
 * be merged to (1,10), but (1,5) and (6,10) would not be merged. (1,5) and (4,10) would also not be merged since they
 * overlap.
 * @param intervals_map
 * @return new ContigToIntervals object
 */
ContigToIntervals MergeAdjacentIntervals(const ContigToIntervals& intervals_map) {
  ContigToIntervals ret;
  for (const auto& [contig, vec] : intervals_map) {
    if (vec.empty()) {
      continue;
    }
    ret[contig].emplace_back(*vec.begin());
    auto ret_it = std::prev(ret[contig].end());
    for (auto it = std::next(vec.begin()); it != vec.end(); ++it) {
      if (it->start == ret_it->end) {
        ret_it->end = it->end;
      } else {
        ret[contig].emplace_back(*it);
        ret_it = std::prev(ret[contig].end());
      }
    }
  }
  return ret;
}

/**
 * @brief given a set of "on-target" intervals, define a set of off-targets by subtracting the on-targets intervals
 * from whole contigs
 * @param intervals_map  - set of on-target intervals
 * @param seq_lens  - map defining lengths of contigs
 * @return
 */
ContigToIntervals GenerateCandidateOffTargets(const ContigToIntervals& intervals_map,
                                              const std::unordered_map<std::string, size_t>& seq_lens) {
  // off-targets = all intervals that exclude on targets
  ContigToIntervals ret;
  // gather each set of intervals under their shared contig
  for (const auto& [contig, length] : seq_lens) {
    auto it = intervals_map.find(contig);
    if (it != intervals_map.end()) {
      std::vector<std::pair<size_t, size_t>> vec;
      size_t p = 0;
      for (const auto& intv : it->second) {
        if (intv.start < p) {
          throw error::Error("Off-target generation failed: overlapping interval on {}:{}-{} (expected start >= {})",
                             contig,
                             intv.start,
                             intv.end,
                             p);
        }
        ret[contig].emplace_back(p, intv.start);
        p = intv.end;
      }
      ret[contig].emplace_back(p, length);
    }
  }
  return ret;
}

/**
 * @brief removes 0-mappability regions from a set of intervals
 * @param map - set of intervals from which to remove 0-mappability intervals
 * @param mapp_bigwig - BigWig object of mappability statistics
 * @return
 */
ContigToIntervals Remove0MappableRegions(const ContigToIntervals& map, BigWig& mapp_bigwig) {
  ContigToIntervals ret;
  for (const auto& [contig, vec] : map) {
    for (auto intv : vec) {
      BwOverlappingIntervalsPtr bw_intvs = mapp_bigwig.GetInterval(contig, intv.start, intv.end);
      if (!bw_intvs || bw_intvs->l == 0) {
        Logging::Warn(
            "could not retrieve region {}:{}-{} from BigWig! skipping this interval...", contig, intv.start, intv.end);
        continue;
      }
      for (size_t i = 0; i < bw_intvs->l; ++i) {
        size_t bw_s = bw_intvs->start[i];
        float val = bw_intvs->value[i];
        if (math::IsCloseToZero(val)) {
          if (intv.start <= bw_s) {
            ret[contig].emplace_back(intv.start, bw_s);
          }
          if (i + 1 < bw_intvs->l) {
            intv.start = bw_intvs->start[i + 1];
          }
        }
      }
      // if the last sub-interval has a value of 0, we're already done
      if (!math::IsCloseToZero(bw_intvs->value[bw_intvs->l - 1])) {
        ret[contig].emplace_back(intv.start, intv.end);
      }
    }
  }
  return ret;
}

/**
 * @brief removes "low" mappability regions from a set of intervals.
 * @param map - intervals from which to remove low mappability
 * @param mapp_bigwig - mappability bigwig object
 * @param min_mappability - minimum threshold of mappability
 * @param sex_chrom_behavior - for setting behavior of how sex chromosomes should be handled - see utility-functions.h
 * @return
 */
ContigToIntervals RemoveLowMappableRegions(const ContigToIntervals& map,
                                           BigWig& mapp_bigwig,
                                           f64 min_mappability,
                                           SexChromHandling sex_chrom_behavior) {
  ContigToIntervals ret;
  for (const auto& [contig, vec] : map) {
    if ((sex_chrom_behavior == SexChromHandling::kNoSexChrom && IsInAllosome(contig)) ||
        (sex_chrom_behavior == SexChromHandling::kOnlySexChrom && !IsInAllosome(contig))) {
      ret[contig] = vec;
      continue;
    }  // else kAllChrom, no skipping
    for (const auto& intv : vec) {
      f64 mean_mappability = mapp_bigwig.GetMean(contig, intv.start, intv.end);
      if (mean_mappability >= min_mappability) {
        ret[contig].emplace_back(intv.start, intv.end);
      }
      // first split the interval by boundaries of 0 mappability regions
    }
  }
  return ret;
}

/**
 * @brief remove intervals not exceeding a given length
 * @param map - intervals
 * @param length - minimum length
 * @return
 */
ContigToIntervals RemoveRegionsByLength(const ContigToIntervals& map, size_t length) {
  ContigToIntervals ret;
  for (const auto& [contig, vec] : map) {
    for (const auto& v : vec) {
      if (v.end - v.start >= length) {
        ret[contig].emplace_back(v);
      }
    }
  }
  return ret;
}

/**
 * @brief split a std::vector of regions into smaller windows with maximum length <w>. If a contig does not divide by
 * <w> perfectly, then <w> is reduced until it is a perfect factor of the contig length.
 * @param intervals_map intervals to split
 * @param w maximum desired width of each window
 * @return a new intervals map
 */
ContigToIntervals SplitIntervalsByWindowSizeAdjustForEqualWindows(const ContigToIntervals& map, size_t interval_size) {
  ContigToIntervals ret;
  for (const auto& [contig, vec] : map) {
    for (const auto& v : vec) {
      size_t length = v.end - v.start;
      // using the width, find how many partitions we should make
      auto n_parts = static_cast<size_t>(std::ceil(static_cast<f64>(length) / static_cast<f64>(interval_size)));
      // determine size of each partition depending on whether the start/end of current window is odd/even
      for (size_t i = 0; i < n_parts; ++i) {
        size_t start = v.start + i * length / n_parts;
        size_t end = v.start + (i + 1) * length / n_parts;
        ret[contig].emplace_back(start, end);
      }
    }
  }
  return ret;
}

/**
 * @brief convert a bed file to a ContigToIntervals map
 * @param bed_file
 * @return
 */
ContigToIntervals BedToMap(const fs::path& bed_file) {
  ContigToIntervals ret;
  std::vector<io::BedRegion> bed_vec(io::ParseBedFile(bed_file));
  for (const auto& bed_region : bed_vec) {
    ret[bed_region.chromosome].emplace_back(bed_region.start, bed_region.end);
  }
  return ret;
}

/**
 * @brief get the number of intervals in a ContigToIntervals object
 * @param map
 * @return
 */
size_t GetNContigToIntervals(const ContigToIntervals& map) {
  size_t n = 0;
  for (const auto& [contig, vec] : map) {
    n += vec.size();
  }
  return n;
}

/**
 * @brief generate a final set of off-targets given a set of on-targets.
 * 1) initial set of candidates is the on-targets subtracted from the whole genome
 * 2) remove 0 mappability regions from candidates
 * 3) remove candidates from small contigs and candidates that cannot be trimmed
 * 4) trim candidates
 * 5) split candidates into windows of specified length
 * 6) remove short candidates
 * 7) remove low mappability candidates
 * @param on_targets - on target intervals
 * @param seq_lengths - object containing lengths of each contig
 * @param mapp_bigwig  - bigwig object containing mappability values for whole genome
 * @param off_target_min_width  - mininum width for off target
 * @param off_target_interval_size  - average window size by which to split off-target candidates
 * @param off_target_trim_length  - length to trim off of each flank of off target candidates
 * @param off_target_min_mappability  - minimum mappability for off target candidates to keep (excluding sex
 * chromosomes)
 * @param sex_target_min_mappability - minimum mappability for any candidate intervals in sex chromosomes
 * @return
 */
ContigToIntervals GenerateOffTargets(const ContigToIntervals& on_targets,
                                     const SeqLenMap& seq_lengths,
                                     BigWig& mapp_bigwig,
                                     const ContigToIntervals& blocklist_intervals,
                                     size_t off_target_min_width,
                                     size_t off_target_interval_size,
                                     size_t off_target_trim_length,
                                     f64 off_target_min_mappability,
                                     f64 sex_target_min_mappability) {
  // remove on target regions
  Logging::Info("generating initial candidates for off-target intervals");
  ContigToIntervals off_targets = GenerateCandidateOffTargets(on_targets, seq_lengths);
  Logging::Info("{} off-target intervals", GetNContigToIntervals(off_targets));
  Logging::Info("Removing blocklisted regions from off-target intervals");
  off_targets = SubtractRegionsFromIntervals(off_targets, blocklist_intervals);
  Logging::Info("{} off-target intervals after removing blocklisted regions", GetNContigToIntervals(off_targets));
  off_targets = Remove0MappableRegions(off_targets, mapp_bigwig);
  // remove entire contigs that are too small
  Logging::Info("{} off-target intervals after removing 0-mappable regions", GetNContigToIntervals(off_targets));
  Logging::Info("removing off-target intervals from small contigs", GetNContigToIntervals(off_targets));
  for (const auto& key_pair : off_targets) {
    auto it = seq_lengths.find(key_pair.first);
    if (it == seq_lengths.end() || it->second < off_target_min_width) {
      off_targets.erase(it->first);
    }
  }
  // remove unmappable regions (mappability == 0)
  // NOTE: this doesn't really happen, so let's just skip this step. Low mappability will be dealt with later.
  // remove regions that cannot be trimmed
  Logging::Info("{} off-target intervals", GetNContigToIntervals(off_targets));
  Logging::Info("removing untrimmable off-target intervals (trim length = {})", off_target_trim_length);
  off_targets = RemoveRegionsByLength(off_targets, 2 * off_target_trim_length);
  // trim off-target regions
  for (auto& [contig, vec] : off_targets) {
    for (auto& v : vec) {
      v.start += off_target_trim_length;
      v.end -= off_target_trim_length;
    }
  }
  // split off-target regions
  Logging::Info("{} off-target intervals", GetNContigToIntervals(off_targets));
  Logging::Info("splitting off-target intervals into windows of {}bp", off_target_interval_size);
  off_targets = SplitIntervalsByWindowSizeAdjustForEqualWindows(off_targets, off_target_interval_size);
  // remove small off-target regions
  Logging::Info("{} off-target intervals", GetNContigToIntervals(off_targets));
  Logging::Info("removing off-target intervals < {}bp", off_target_min_width);
  off_targets = RemoveRegionsByLength(off_targets, off_target_min_width);
  // remove low mappability regions
  Logging::Info("{} off-target intervals", GetNContigToIntervals(off_targets));
  Logging::Info("removing off-target intervals with <{} mappability", off_target_min_mappability);
  off_targets =
      RemoveLowMappableRegions(off_targets, mapp_bigwig, off_target_min_mappability, SexChromHandling::kNoSexChrom);
  off_targets =
      RemoveLowMappableRegions(off_targets, mapp_bigwig, sex_target_min_mappability, SexChromHandling::kOnlySexChrom);
  Logging::Info("{} off-target intervals", GetNContigToIntervals(off_targets));
  // merge in off-target regions
  return off_targets;
}

/**
 * @brief converts ContigToIntervals to a std::vector of region strings
 * @param map
 * @return
 */
std::vector<std::string> ToRegionStringVec(const ContigToIntervals& map) {
  std::vector<std::string> ret;
  for (const auto& [contig, vec] : map) {
    for (const auto& v : vec) {
      ret.emplace_back(fmt::format("{}:{}-{}", contig, v.start + 1, v.end));
    }
  }
  return ret;
}

/**
 * @brief given a map of contig->interval-lists, generate evenly sized intervals
 * across the sequences. If <interval_size> does not divide perfectly into the
 * size of a contig, then th size of the last interval over that contig will be
 * the size of the remainder after dividing the contig length by the interval
 * size
 * @param seq_lengths SeqLenMap of sequences/lengths from which to derive intervals
 * @param interval_size  desired interval size
 * @return
 */
ContigToIntervals SplitIntervalsByWindowSizedLastWindowSizeOfRemainder(const ContigToIntervals& map,
                                                                       size_t interval_size) {
  ContigToIntervals ret;
  for (const auto& [contig, vec] : map) {
    for (const auto& v : vec) {
      for (size_t i = v.start; i < v.end; i += interval_size) {
        size_t end = i + interval_size > v.end ? v.end : i + interval_size;
        ret[contig].emplace_back(i, end);
      }
    }
  }
  return ret;
}

/**
 * @brief given a map of contig->lengths, generate intervals that span the entire length of each contig
 * @param seq_lengths SeqLenMap of sequences/lengths from which to derive intervals
 * @param include_alt_contigs Boolean to indicate whether to include non-canonical contigs
 */
ContigToIntervals IntervalsFromSeqLenMap(const SeqLenMap& seq_lengths, bool include_alt_contigs) {
  ContigToIntervals ret;
  std::regex pattern("(_|chrEBV|chrM)");
  for (const auto& [contig, len] : seq_lengths) {
    if (!include_alt_contigs && std::regex_search(contig, pattern)) {
      Logging::Info("Skipping non-canonical contig: {}", contig);
      continue;
    }
    ret[contig].emplace_back(0, len);
  }
  return ret;
}

void ThrowIfEmpty(const ContigToIntervals& map, const std::optional<std::string>& diagnostic_message) {
  if (map.empty() || std::all_of(map.begin(), map.end(), [](const auto& v) { return v.second.empty(); })) {
    if (diagnostic_message) {
      throw std::runtime_error(*diagnostic_message);
    } else {
      throw std::runtime_error("empty ContigToIntervalsMap");
    }
  }
}

std::vector<Interval> SubtractIntervalBFromIntervalA(const Interval& a, const Interval& b) {
  if (a.start >= b.end || a.end <= b.start) {
    // no overlap, keep a in its entirety
    return {a};
  }
  if (b.start <= a.start && b.end >= a.end) {
    // b contains a
    return {};
  }
  if (a.start <= b.start && a.end >= b.end) {
    // a contains b
    return {{a.start, b.start}, {b.end, a.end}};
  }
  if (a.start < b.start) {
    // a overlaps left side of b
    return {{a.start, b.start}};
  }
  if (a.end > b.end) {
    // a overlaps right side of b
    return {{b.end, a.end}};
  }
  //  other cases should not happen
  Logging::Warn("interval subtraction {}-{} from {}-{} is not handled", b.start, b.end, a.start, a.end);
  return {};
}

/**
 * @brief subtract intervals in `intervals_to_remove` from `intervals`. Assumes that intervals in both sets do not
 * overlap with one another.
 * @param intervals
 * @param intervals_to_remove
 * @return a new set of intervals representing the original `intervals` set "minus" the `intervals_to_remove` set
 */
ContigToIntervals SubtractRegionsFromIntervals(const ContigToIntervals& intervals,
                                               const ContigToIntervals& intervals_to_remove) {
  ContigToIntervals ret;
  for (const auto& [contig, intvs] : intervals) {
    ret[contig] = {};
    std::vector<Interval> intvs_to_remove;
    if (intervals_to_remove.find(contig) == intervals_to_remove.end()) {
      // no intervals to remove, just copy the current intervals
      ret[contig] = intvs;
      continue;
    } else {
      intvs_to_remove = intervals_to_remove.at(contig);
    }
    auto start_it = intvs_to_remove.begin();
    for (const auto& intv : intvs) {
      // first, disregard any interval in intervals_to_remove that ends before the current interval starts
      // the logic here also lets us keep intvs in intervals_to_remove that might overlap 2+ intvs in intervals
      while (start_it != intvs_to_remove.end() && start_it->end < intv.start) {
        ++start_it;
      }
      if (start_it == intvs_to_remove.end()) {
        // no more intervals to remove, just copy the current interval
        ret[contig].emplace_back(intv);
        continue;
      }
      // then subtract any interval in intervals_to_remove that starts before the current interval ends
      // assume that none of the intervals in intervals_to_remove overlap with one another
      std::deque<Interval> remaining_intervals{intv};
      for (auto it = start_it; it != intvs_to_remove.end(); ++it) {
        // if the remaining_interval occurs before *it, then move it to the result
        while (!remaining_intervals.empty() && remaining_intervals.front().end < it->start) {
          ret[contig].emplace_back(remaining_intervals.front());
          remaining_intervals.pop_front();
        }
        // subtract *it from the remaining interval
        // this should always return 1 or 2 intervals
        if (!remaining_intervals.empty()) {
          auto subtracted_intvs = SubtractIntervalBFromIntervalA(remaining_intervals.front(), *it);
          remaining_intervals.pop_front();
          for (const auto& sub_intv : subtracted_intvs) {
            remaining_intervals.push_back(sub_intv);
          }
        }
      }
      // push any remaining intervals that were not subtracted
      while (!remaining_intervals.empty()) {
        ret[contig].emplace_back(remaining_intervals.front());
        remaining_intervals.pop_front();
      }
    }
  }
  return ret;
}

}  // namespace xoos::cnc
