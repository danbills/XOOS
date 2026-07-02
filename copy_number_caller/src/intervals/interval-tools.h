#pragma once
#include <optional>
#include <string>
#include <vector>

#include <xoos/io/bed-region.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>

#include "bigwig_cpp.h"
#include "intervals/intervals.h"
#include "io/fai.h"
#include "sex.h"

namespace xoos::cnc {

ContigToIntervals MergeAdjacentIntervals(const ContigToIntervals& intervals_map);
ContigToIntervals GenerateCandidateOffTargets(const ContigToIntervals& intervals_map,
                                              const std::unordered_map<std::string, size_t>& seq_lens);
ContigToIntervals RemoveLowMappableRegions(const ContigToIntervals& map,
                                           BigWig& mapp_bigwig,
                                           f64 min_mappability,
                                           SexChromHandling sex_chrom_behavior = SexChromHandling::kNoSexChrom);
ContigToIntervals RemoveRegionsByLength(const ContigToIntervals& map, size_t length);
ContigToIntervals SplitIntervalsByWindowSizeAdjustForEqualWindows(const ContigToIntervals& map, size_t interval_size);
ContigToIntervals BedToMap(const fs::path& bed_file);
ContigToIntervals GenerateOffTargets(const ContigToIntervals& on_targets,
                                     const SeqLenMap& seq_lengths,
                                     BigWig& mapp_bigwig,
                                     const ContigToIntervals& blocklist_intervals,
                                     size_t off_target_min_width,
                                     size_t off_target_interval_size,
                                     size_t off_target_trim_length,
                                     f64 off_target_min_mappability,
                                     f64 sex_target_min_mappability);
ContigToIntervals MergeIntervals(const ContigToIntervals& map1, const ContigToIntervals& map2);
size_t GetNContigToIntervals(const ContigToIntervals& map);
std::vector<std::string> ToRegionStringVec(const ContigToIntervals& map);
ContigToIntervals Remove0MappableRegions(const ContigToIntervals& map, BigWig& mapp_bigwig);
ContigToIntervals SplitIntervalsByWindowSizedLastWindowSizeOfRemainder(const ContigToIntervals& map,
                                                                       size_t interval_size);
ContigToIntervals IntervalsFromSeqLenMap(const SeqLenMap& seq_lengths, bool include_alt_contigs);
std::vector<Interval> SubtractIntervalBFromIntervalA(const Interval& a, const Interval& b);
ContigToIntervals SubtractRegionsFromIntervals(const ContigToIntervals& intervals,
                                               const ContigToIntervals& intervals_to_remove);
void ThrowIfEmpty(const ContigToIntervals& map, const std::optional<std::string>& diagnostic_message = std::nullopt);
}  // namespace xoos::cnc
