#include "parallel-compute-utils.h"

#include <algorithm>
#include <limits>
#include <mutex>
#include <string_view>
#include <unordered_set>

#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/for_each.hpp>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/util/math.h>

#include "compute-bam-features/region.h"
#include "compute-vcf-features/compute-vcf-features.h"
#include "compute-vcf-features/vcf-header-util.h"
#include "core/vcf-fields.h"
#include "vcf-to-bed/vcf-to-bed.h"

namespace xoos::svc {

WorkerContext::WorkerContext(const std::string& vcf_file,
                             const std::optional<std::string>& popaf_file,
                             const std::vector<fs::path>& bam_inputs,
                             AlignmentReaderCache& alignment_reader_cache)
    // Set up pointers to BAM readers, SAM headers, and BAM indexes for all input BAM files
    : filter_variants_reader{vcf_file}, vcf_feature_extraction_reader{vcf_file} {
  alignment_readers = alignment_reader_cache.Open(bam_inputs, "r");
  // Set up the population allele frequency VCF reader if provided
  if (popaf_file.has_value()) {
    popaf_reader = io::VcfReader(popaf_file.value());
  }
}

/**
 * @brief Partition a VCF file into regions based on the number of threads available.
 * Each region will contain approximately the same number of variants, allowing for parallel processing.
 * @param vcf_file Path to the VCF file to partition.
 * @param threads Number of threads to use for parallel processing.
 * @return A vector of TargetRegion objects representing the partitioned regions.
 */
vec<TargetRegion> PartitionVcfRegions(const fs::path& vcf_file,
                                      const size_t threads,
                                      std::optional<ChromIntervalsMap> bed_regions) {
  vec<TargetRegion> partitioned_regions;
  const auto safe_threads = std::max<size_t>(threads, 1);
  const u32 left_pad = 1;
  const u32 right_pad = 1;
  ChromIntervalsMap target_regions;
  if (bed_regions.has_value()) {
    target_regions = bed_regions.value();
  }

  // Only collapse intervals if they are immediately adjacent to each other. So, each chromosome would have
  // approximately one interval for each variant. The only exception is large deletions, which have large intervals.
  const u32 collapse_dist = 0;

  // Use a single thread to extract variant positions. If the VCF file has many chromosomes, then the parallelized
  // version (1 parallel task per chromomosome) may cause this step to become very slow.
  auto chrom_to_intervals =
      ExtractVariantIntervalsSingleThreaded(vcf_file, target_regions, left_pad, right_pad, collapse_dist);
  // The goal here is to assign roughly same number of variants per parallel task in each chromosome. Since features
  // have already been computed, there is no need to create parallel tasks for small equal-sized regions across the
  // entire genome. Otherwise, many tasks would be created for regions that do not overlap any variants at all.
  for (auto& [chrom, intervals] : chrom_to_intervals) {
    if (!intervals.empty()) {
      auto num_intervals = intervals.size();
      if (num_intervals > safe_threads) {
        // All partitions of this chromosome would contain near-equal number of intervals.
        const auto num_intervals_per_thread = (num_intervals + safe_threads - 1) / safe_threads;
        auto last_i = num_intervals - 1;
        for (size_t i = 0; i < num_intervals; i += num_intervals_per_thread) {
          auto max_i = std::min(last_i, i + num_intervals_per_thread - 1);
          partitioned_regions.emplace_back(chrom, intervals[i].start, intervals[max_i].end);
        }
      } else {
        // Too few intervals in this chromosome.
        // Create only one partition for all intervals in this chromosome.
        partitioned_regions.emplace_back(chrom, intervals[0].start, (intervals.end() - 1)->end);
      }
    }
  }
  return partitioned_regions;
}

/**
 * @brief Computes BAM and VCF features for a given region.
 * @param global_ctx Global context containing configuration and reference sequences.
 * @param worker_ctx Worker context containing BAM readers, SAM headers, and BAM indexes.
 * @param region The target region for which to compute features.
 * @param bed_regions Optional map of chromosome to vector of intervals representing BED regions.
 * @param interest_regions Optional map of chromosome to vector of intervals representing regions of interest.
 * @return A tuple containing the computed VCF features and BAM features.
 */
std::tuple<VarIdToVcfFeatures, BamRegionFeatureCollection> ComputeBamAndVcfFeaturesForRegion(
    const GlobalContext& global_ctx,
    WorkerContext& worker_ctx,
    const TargetRegion& region,
    const ChromIntervalsMap& bed_regions,
    const ChromIntervalsMap& interest_regions) {
  // This function is designed to be used by an individual worker thread as part of a larger taskflow workflow where
  // each task is getting features for a given region.

  // Compute VCF features first. Once VCF features have been computed we can use the feature set when computing BAM
  // features. This should reduce peak memory usage and lower compute time.
  // If a set of BED Regions has been passed use the intersection of the target regions and the BAM region
  auto bam_region = Region{.chrom = region.chrom, .start = region.start, .end = region.end};
  ChromIntervalsMap target_regions;
  auto i_region = Interval{.start = region.start, .end = region.end};
  if (bed_regions.contains(region.chrom)) {
    for (const auto& interval : bed_regions.at(region.chrom)) {
      if (IntervalOverlap(interval, i_region)) {
        // Adjust interval if it falls partly outside the region
        auto interval_to_use = Interval{.start = interval.start, .end = interval.end};
        if (interval_to_use.start < i_region.start) {
          interval_to_use.start = i_region.start;
        }
        if (interval_to_use.end > i_region.end) {
          interval_to_use.end = i_region.end;
        }
        target_regions[region.chrom].emplace_back(interval_to_use);
      }
    }
  } else if (bed_regions.empty()) {
    target_regions[region.chrom].emplace_back(i_region);
  }
  VarIdToVcfFeatures vcf_features;
  BamRegionFeatureCollection bam_features;
  ChromPosToRefBamFeatures ref_features;
  bool use_vcf_features = true;
  if (!global_ctx.model_config.HasVcfFeatureScoringCols() || !global_ctx.model_config.use_vcf_features) {
    use_vcf_features = false;
  }
  if (!target_regions.empty()) {
    std::optional<StrUnorderedSet> skip_variants = std::nullopt;
    if (global_ctx.skip_variants_vcf.has_value()) {
      skip_variants = std::make_optional(ExtractVariantKeySet(global_ctx.skip_variants_vcf.value(), bam_region));
    }

    // Only compute VCF features for the region if there is a target entry and we have to based on workflow and scoring
    // cols
    std::optional<VarIdToVcfFeatures> passed_vcf_features = std::nullopt;
    if (use_vcf_features) {
      const VcfFeatureExtractionParams extraction_params{
          global_ctx.is_germline_tagging, global_ctx.pre_scan_variants, global_ctx.pre_scan_density};
      vcf_features = ExtractFeaturesForRegion(worker_ctx.vcf_feature_extraction_reader,
                                              global_ctx.genome,
                                              worker_ctx.popaf_reader,
                                              target_regions,
                                              interest_regions,
                                              region.chrom,
                                              extraction_params);
      passed_vcf_features = vcf_features;
    }

    if (global_ctx.ref_seqs == nullptr) {
      throw error::Error("GlobalContext::ref_seqs is null — SetReferenceSequences() was not called");
    }
    bam_features =
        ComputeBamRegionFeatures(worker_ctx.alignment_readers, global_ctx.bam_feat_params, nullptr)
            .ComputeBamFeatures(bam_region, global_ctx.ref_seqs->at(region.chrom), passed_vcf_features, skip_variants);
  }
  return std::make_tuple(vcf_features, bam_features);
}

/**
 * @brief Defines a region query for a single chromosome during the parallel VCF pre-scan.
 */
struct PreScanChromTask {
  // Chromosome name
  std::string chrom;
  // htslib contig index (tid)
  s32 chrom_index;
  // 0-based start position of the region to scan
  u64 start;
  // End position (exclusive) of the region to scan
  u64 end;
};

/**
 * @brief Identifies which VCF samples to extract DP from and how they map to normal/tumor roles.
 * @details sample_indexes holds the 0-based VCF sample indexes to query. normal_pos and tumor_pos
 * are positions within sample_indexes that correspond to the normal and tumor samples respectively
 * (-1 when the role is absent). This allows ExtractDpValues to work with any number of samples
 * while still supporting the downstream ChromMedianDepth normal/tumor split.
 */
struct SampleIndexInfo {
  // 0-based VCF sample indexes to extract DP from
  vec<s32> sample_indexes;
  // Position within sample_indexes for the normal sample (-1 if absent)
  s32 normal_pos{-1};
  // Position within sample_indexes for the tumor sample (-1 if absent)
  s32 tumor_pos{-1};
};

vec<u32> ComputeVariantDensity(const vec<VariantPreScanEntry>& variants) {
  const auto num_variants = variants.size();
  // Each variant counts itself
  vec<u32> density(num_variants, 1);
  for (u64 i = 0; i < num_variants; ++i) {
    const auto max_pos = variants[i].pos + kPreScanDensityWindow;
    for (u64 j = i + 1; j < num_variants && variants[j].pos <= max_pos; ++j) {
      ++density[i];
      ++density[j];
    }
  }
  return density;
}

/**
 * @brief Extract DP values from a VCF record for each sample in the index list.
 * @param record VCF record to extract DP from
 * @param sample_info Sample index information specifying which VCF samples to query
 * @return Vector of per-sample DP values, one per entry in sample_info.sample_indexes
 */
static vec<u32> ExtractDpValues(const io::VcfRecordPtr& record, const SampleIndexInfo& sample_info) {
  const auto& dp_values = record->GetFormatFieldNoCheck<s32>(kFieldDp);
  vec<u32> dps;
  dps.reserve(sample_info.sample_indexes.size());
  for (const auto idx : sample_info.sample_indexes) {
    if (0 <= idx && std::cmp_less(idx, dp_values.size())) {
      dps.emplace_back(static_cast<u32>(dp_values.at(static_cast<size_t>(idx))));
    } else {
      dps.emplace_back(0);
    }
  }
  return dps;
}

/**
 * @brief Check whether a VCF record at the given position passes BED interval filtering.
 * Advances the BED iterator past intervals that end before the current position.
 * @return true if the record passes (should be included), false if it should be skipped
 */
static bool PassesBedFilter(const u64 pos,
                            const std::string_view ref_allele,
                            vec<Interval>& target_intervals,
                            vec<Interval>::iterator& target_itr) {
  if (target_intervals.empty() || target_itr == target_intervals.end()) {
    return false;
  }
  while (target_itr != target_intervals.end() && target_itr->end < pos) {
    ++target_itr;
  }
  if (target_itr == target_intervals.end()) {
    return false;
  }
  return IntervalOverlap(*target_itr, pos, pos + ref_allele.length());
}

/**
 * @brief Build per-chromosome pre-scan tasks from VCF contigs and optional BED regions.
 * @details When BED regions are provided, each task spans the BED extent for that chromosome.
 * Otherwise, each task spans the full contig length (or s64 max if the length is missing).
 * @param vcf_reader Open VCF reader used to query contig indexes and lengths
 * @param bed_regions Optional BED intervals to restrict the scan range per chromosome
 * @return Vector of PreScanChromTask, one per chromosome to scan
 */
static vec<PreScanChromTask> BuildPreScanTasks(io::VcfReader& vcf_reader,
                                               const std::optional<ChromIntervalsMap>& bed_regions) {
  vec<PreScanChromTask> tasks;
  const auto hdr = vcf_reader.GetHeader();
  const auto contig_indexes = vcf_reader.GetContigIndexes();
  const auto contig_lengths = hdr->GetContigLengths();

  if (bed_regions.has_value() && !bed_regions->empty()) {
    for (const auto& [chrom, intervals] : bed_regions.value()) {
      if (contig_indexes.contains(chrom) && !intervals.empty()) {
        const auto cid = contig_indexes.at(chrom);
        tasks.emplace_back(chrom, cid, intervals.begin()->start, (intervals.end() - 1)->end);
      }
    }
  } else {
    for (const auto& [chrom, cid] : contig_indexes) {
      // Use contig length from the header when available; fall back to signed-int max so
      // htslib scans the entire contig even when ##contig length metadata is missing or zero.
      // s64 max is used because VcfReader::SetRegion forwards to htslib's hts_pos_t (int64_t).
      const u64 end = contig_lengths.contains(chrom) && contig_lengths.at(chrom) > 0
                          ? contig_lengths.at(chrom)
                          : static_cast<u64>(std::numeric_limits<s64>::max());
      tasks.emplace_back(chrom, cid, 0, end);
    }
  }
  return tasks;
}

/**
 * @brief Scan a single chromosome region for variant DP and position data using an indexed VCF query.
 * @details Opens its own VcfReader so it can be called from a parallel thread without sharing state.
 * De-duplicates by position and optionally filters by BED intervals.
 * @param vcf_file Path to the input VCF file
 * @param task Chromosome region to scan (chrom, index, start, end)
 * @param sample_info Sample index information for DP extraction
 * @param has_bed Whether BED filtering is active
 * @param bed BED intervals keyed by chromosome (only consulted when has_bed is true)
 * @return Sorted vector of VariantPreScanEntry for the scanned region
 */
static vec<VariantPreScanEntry> ScanChromosomeVariants(const fs::path& vcf_file,
                                                       const PreScanChromTask& task,
                                                       const SampleIndexInfo& sample_info,
                                                       const bool has_bed,
                                                       const ChromIntervalsMap& bed) {
  vec<VariantPreScanEntry> entries;
  std::unordered_set<u64> seen_positions;

  try {
    io::VcfReader vcf_reader(vcf_file);
    if (!vcf_reader.SetRegion(task.chrom_index, task.start, task.end)) {
      return entries;
    }

    vec<Interval> target_intervals;
    vec<Interval>::iterator target_itr;
    const bool has_chrom_bed = has_bed && bed.contains(task.chrom);
    if (has_chrom_bed) {
      target_intervals = bed.at(task.chrom);
      target_itr = target_intervals.begin();
    }

    while (const auto& record = vcf_reader.GetNextRegionRecord(BCF_UN_ALL)) {
      if (record->Position() < 0) {
        continue;
      }
      const auto pos = static_cast<u64>(record->Position());

      if (has_chrom_bed && !PassesBedFilter(pos, record->Allele(0), target_intervals, target_itr)) {
        continue;
      }

      if (!seen_positions.insert(pos).second) {
        continue;
      }

      const auto ref_len = static_cast<u32>(record->Allele(0).length());
      auto dps = ExtractDpValues(record, sample_info);
      u32 total = 0;
      for (const auto dp : dps) {
        total += dp;
      }
      entries.emplace_back(pos, total, std::move(dps), ref_len);
    }
  } catch (const std::exception& e) {
    throw error::Error("Error in VCF pre-scan for chromosome '{}': {}", task.chrom, e.what());
  }

  std::ranges::sort(entries, [](const VariantPreScanEntry& a, const VariantPreScanEntry& b) { return a.pos < b.pos; });
  return entries;
}

/**
 * @brief Single-threaded sequential VCF scan fallback used when no index is available or only one thread is requested.
 * @details Iterates all records in file order, grouping entries by chromosome. De-duplicates by position and
 * optionally filters by BED intervals.
 * @param vcf_file Path to the input VCF file
 * @param sample_info Sample index information for DP extraction
 * @param has_bed Whether BED filtering is active
 * @param bed BED intervals keyed by chromosome (only consulted when has_bed is true)
 * @return Map of chromosome name to sorted vector of VariantPreScanEntry
 */
static StrMap<vec<VariantPreScanEntry>> ScanVariantsSequential(const fs::path& vcf_file,
                                                               const SampleIndexInfo& sample_info,
                                                               const bool has_bed,
                                                               const ChromIntervalsMap& bed) {
  StrMap<vec<VariantPreScanEntry>> chrom_variants;
  const io::VcfReader vcf_reader(vcf_file);
  std::string prev_chrom;
  vec<VariantPreScanEntry> entries;
  std::unordered_set<u64> seen_positions;

  vec<Interval> target_intervals;
  vec<Interval>::iterator target_itr;

  while (const auto& record = vcf_reader.GetNextRecord()) {
    if (record->Position() < 0) {
      continue;
    }
    const std::string& chrom = record->Chromosome();
    const auto pos = static_cast<u64>(record->Position());

    if (prev_chrom != chrom) {
      if (!entries.empty()) {
        chrom_variants[prev_chrom] = std::move(entries);
        entries = {};
      }
      seen_positions.clear();
      prev_chrom = chrom;
      if (has_bed && bed.contains(chrom)) {
        target_intervals = bed.at(chrom);
        target_itr = target_intervals.begin();
      } else {
        target_intervals = {};
      }
    }

    if (has_bed && !PassesBedFilter(pos, record->Allele(0), target_intervals, target_itr)) {
      continue;
    }

    if (!seen_positions.insert(pos).second) {
      continue;
    }

    const auto ref_len = static_cast<u32>(record->Allele(0).length());
    auto dps = ExtractDpValues(record, sample_info);
    u32 total = 0;
    for (const auto dp : dps) {
      total += dp;
    }
    entries.emplace_back(pos, total, std::move(dps), ref_len);
  }
  if (!entries.empty()) {
    chrom_variants[prev_chrom] = std::move(entries);
  }
  return chrom_variants;
}

/**
 * @brief Compute per-chromosome median depths from collected variant DP values.
 * @param chrom_variants Map of chromosome name to pre-scan entries containing per-variant DP values
 * @param sample_info Sample index information with normal/tumor position mappings
 * @return ChromMedianDepth with normal and tumor median DP per chromosome
 */
static ChromMedianDepth ComputeChromMedianDepths(const StrMap<vec<VariantPreScanEntry>>& chrom_variants,
                                                 const SampleIndexInfo& sample_info) {
  ChromMedianDepth chrom_median_dp;
  for (const auto& [chrom, entries] : chrom_variants) {
    vec<u32> normal_dps;
    vec<u32> tumor_dps;
    normal_dps.reserve(entries.size());
    tumor_dps.reserve(entries.size());
    for (const auto& entry : entries) {
      if (sample_info.normal_pos >= 0 && std::cmp_less(sample_info.normal_pos, entry.sample_dps.size())) {
        normal_dps.emplace_back(entry.sample_dps.at(static_cast<size_t>(sample_info.normal_pos)));
      }
      if (sample_info.tumor_pos >= 0 && std::cmp_less(sample_info.tumor_pos, entry.sample_dps.size())) {
        tumor_dps.emplace_back(entry.sample_dps.at(static_cast<size_t>(sample_info.tumor_pos)));
      }
    }
    if (!normal_dps.empty()) {
      chrom_median_dp.normal[chrom] = math::Median(normal_dps);
    }
    if (!tumor_dps.empty()) {
      chrom_median_dp.tumor[chrom] = math::Median(tumor_dps);
    }
  }
  return chrom_median_dp;
}

VcfPreScanResult ParallelVcfPreScan(const fs::path& vcf_file,
                                    const size_t threads,
                                    std::optional<ChromIntervalsMap> bed_regions) {
  Logging::Info("Starting parallel VCF pre-scan of {}", vcf_file);

  // Determine sample indexes
  SampleIndexInfo sample_info;
  bool use_single_thread = false;
  vec<PreScanChromTask> tasks;
  {
    io::VcfReader vcf_reader(vcf_file);
    const auto hdr = vcf_reader.GetHeader();

    if (!hdr->HasInfoField(kFieldDp)) {
      Logging::Warn("VCF file does not have FORMAT field '{}'; pre-scan DP values will be zero", kFieldDp);
    }

    const auto tn_sample_idx = GetTumorNormalSampleIndexes(hdr);
    if (tn_sample_idx.has_value()) {
      // Tumor-normal VCF: extract DP from both samples
      sample_info.sample_indexes = {tn_sample_idx->normal_sample_idx, tn_sample_idx->tumor_sample_idx};
      sample_info.normal_pos = 0;
      sample_info.tumor_pos = 1;
    } else {
      // Single-sample VCF: extract DP from sample 0 as normal
      sample_info.sample_indexes = {0};
      sample_info.normal_pos = 0;
    }

    use_single_thread = threads <= 1 || !vcf_reader.HasIndex();
    if (!use_single_thread) {
      tasks = BuildPreScanTasks(vcf_reader, bed_regions);
    }
  }

  const bool has_bed = bed_regions.has_value() && !bed_regions->empty();
  const ChromIntervalsMap empty_bed;
  const auto& bed = has_bed ? bed_regions.value() : empty_bed;

  StrMap<vec<VariantPreScanEntry>> chrom_variants;

  if (use_single_thread) {
    chrom_variants = ScanVariantsSequential(vcf_file, sample_info, has_bed, bed);
  } else {
    Logging::Info("Pre-scanning VCF with {} threads across {} chromosomes", threads, tasks.size());
    std::mutex result_mutex;
    tf::Executor executor(threads);
    tf::Taskflow flow;
    flow.for_each(
            tasks.begin(),
            tasks.end(),
            [&vcf_file, &sample_info, has_bed, &bed, &result_mutex, &chrom_variants](const PreScanChromTask& task) {
              auto entries = ScanChromosomeVariants(vcf_file, task, sample_info, has_bed, bed);
              if (!entries.empty()) {
                const std::scoped_lock lock{result_mutex};
                chrom_variants[task.chrom] = std::move(entries);
              }
            })
        .name("Parallel VCF pre-scan");
    executor.run(flow).get();
  }

  auto chrom_median_dp = ComputeChromMedianDepths(chrom_variants, sample_info);

  StrMap<vec<u32>> chrom_variant_density;
  for (const auto& [chrom, entries] : chrom_variants) {
    chrom_variant_density[chrom] = ComputeVariantDensity(entries);
  }

  Logging::Info("Pre-scan complete: {} chromosomes with variants", chrom_variants.size());

  return VcfPreScanResult{
      .chrom_variants = std::move(chrom_variants),
      .chrom_median_dp = std::move(chrom_median_dp),
      .chrom_variant_density = std::move(chrom_variant_density),
  };
}

/**
 * @brief Filter pre-scan entries to only those overlapping BED intervals for a single chromosome.
 * @details Uses an advancing-iterator pattern over sorted intervals and sorted entries for efficiency.
 * @param entries Pre-scan entries for one chromosome (sorted by position)
 * @param intervals Sorted BED intervals for the same chromosome
 * @return Entries that overlap at least one BED interval
 */
static vec<VariantPreScanEntry> FilterEntriesByBed(const vec<VariantPreScanEntry>& entries,
                                                   const vec<Interval>& intervals) {
  vec<VariantPreScanEntry> filtered;
  auto itr = intervals.begin();
  for (const auto& entry : entries) {
    const u64 var_end = entry.pos + entry.ref_len;
    while (itr != intervals.end() && itr->end < entry.pos) {
      ++itr;
    }
    if (itr != intervals.end() && IntervalOverlap(*itr, entry.pos, var_end)) {
      filtered.emplace_back(entry);
    }
  }
  return filtered;
}

/**
 * @brief Partition a single chromosome's entries into DP-weighted regions.
 * @param chrom Chromosome name
 * @param entries Pre-scan entries for this chromosome (non-empty)
 * @param threads Number of target regions
 * @param left_pad Left padding for the first region boundary
 * @param right_pad Right padding for the last region boundary
 * @param out Output vector to append regions to
 */
static void PartitionChromByDp(const std::string& chrom,
                               const vec<VariantPreScanEntry>& entries,
                               const size_t threads,
                               const u32 left_pad,
                               const u32 right_pad,
                               vec<TargetRegion>& out) {
  if (entries.size() <= threads) {
    const u64 start = entries.front().pos >= left_pad ? entries.front().pos - left_pad : 0;
    const u64 end = entries.back().pos + entries.back().ref_len + right_pad;
    out.emplace_back(chrom, start, end);
    return;
  }

  u64 total_weight = 0;
  for (const auto& entry : entries) {
    total_weight += std::max(entry.total_dp, static_cast<u32>(1));
  }

  const u64 min_target = std::max(total_weight / threads, static_cast<u64>(1));
  u64 accumulated_weight = 0;
  size_t region_start_idx = 0;
  // Maximum reference end (pos + ref_len) of any entry in the current region block. A region
  // boundary is placed at the next entry's position, so cutting while the next entry starts
  // before this end would split a deletion's reference span across two regions and orphan the
  // downstream wildcard-ALT (*) records that depend on it. Track the furthest span reach so
  // the cut can be deferred until it is safe.
  u64 region_max_ref_end = 0;

  for (size_t i = 0; i < entries.size(); ++i) {
    accumulated_weight += std::max(entries[i].total_dp, static_cast<u32>(1));
    region_max_ref_end = std::max(region_max_ref_end, entries[i].pos + entries[i].ref_len);

    const bool is_last = (i == entries.size() - 1);
    // Defer an intermediate cut while the next entry falls strictly inside the current block's
    // deletion span, so deletion-spanned downstream variants stay in the same region.
    const bool span_safe = is_last || entries[i + 1].pos >= region_max_ref_end;
    if ((accumulated_weight >= min_target && span_safe) || is_last) {
      const u64 start = region_start_idx == 0 && entries[region_start_idx].pos >= left_pad
                            ? entries[region_start_idx].pos - left_pad
                            : entries[region_start_idx].pos;

      // For the last region, extend past the final variant. For intermediate regions,
      // set end = next region's first variant position to form non-overlapping half-open
      // intervals [start, end).
      const u64 end = is_last ? region_max_ref_end + right_pad : entries[i + 1].pos;

      out.emplace_back(chrom, start, end);
      region_start_idx = i + 1;
      accumulated_weight = 0;
      region_max_ref_end = 0;
    }
  }
}

vec<TargetRegion> PartitionVcfRegionsByDp(const VcfPreScanResult& pre_scan,
                                          const size_t threads,
                                          const u32 left_pad,
                                          const u32 right_pad,
                                          const std::optional<ChromIntervalsMap>& bed_regions) {
  vec<TargetRegion> partitioned_regions;
  const bool has_bed = bed_regions.has_value() && !bed_regions->empty();

  for (const auto& [chrom, all_entries] : pre_scan.chrom_variants) {
    if (all_entries.empty()) {
      continue;
    }

    if (has_bed) {
      if (!bed_regions->contains(chrom)) {
        continue;
      }
      const auto filtered = FilterEntriesByBed(all_entries, bed_regions->at(chrom));
      if (!filtered.empty()) {
        PartitionChromByDp(chrom, filtered, threads, left_pad, right_pad, partitioned_regions);
      }
    } else {
      PartitionChromByDp(chrom, all_entries, threads, left_pad, right_pad, partitioned_regions);
    }
  }

  return partitioned_regions;
}

}  // namespace xoos::svc
