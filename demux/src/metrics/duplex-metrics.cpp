#include "metrics/duplex-metrics.h"

#include <fmt/format.h>
#include <xoos/histogram/histogram-summary.h>
#include <xoos/log/logging.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <functional>
#include <numeric>
#include <string>
#include <type_traits>

#include "core/demux-and-trim-pipeline.h"
#include "csv.hpp"
#include "metric-filenames.h"
#include "metrics-constants.h"
#include "metrics-constraints.h"
#include "metrics.h"

namespace xoos::demux {

thread_local concurrent::EnumerableThreadLocal<DuplexMetrics> DuplexMetrics::instance{
    std::make_shared<DuplexMetrics>()};

DuplexMetrics::DuplexMetrics()
    : total_length_distr(metrics_constraints::max_logged_read_length + 1, 0),
      unassigned_length_distr(metrics_constraints::max_logged_read_length + 1, 0),
      no_hairpin_length_distr(metrics_constraints::max_logged_read_length + 1, 0),
      full_duplex_length_distr(metrics_constraints::max_sid_id_index + 1,
                               LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)),
      partial_duplex_length_distr(metrics_constraints::max_sid_id_index + 1,
                                  LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)),
      endadapter_position_distr(metrics_constraints::max_sid_id_index + 1,
                                LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)),
      passing_length_distr(metrics_constraints::max_sid_id_index + 1,
                           LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)) {}

static void AddOther(PerSIDHistogram& hist, const PerSIDHistogram& other) {
  for (size_t i = 0; i < hist.size(); ++i) {
    hist[i] += other[i];
  }
}

static void AddOther(PerSIDCount& count, const PerSIDCount& other) {
  for (size_t i = 0; i < count.size(); ++i) {
    count[i] += other[i];
  }
}

DuplexMetrics& DuplexMetrics::Instance() { return *instance.Local(); }

void DuplexMetrics::Add(const DuplexMetrics& other) {
  // merge unassigned counts
  unassigned_counts.read_too_long += other.unassigned_counts.read_too_long;
  unassigned_counts.no_hairpin_found += other.unassigned_counts.no_hairpin_found;
  unassigned_counts.read_too_short += other.unassigned_counts.read_too_short;
  unassigned_counts.raw_bases += other.unassigned_counts.raw_bases;

  // merge failed assigned counts
  AddOther(failed_assigned_counts.consensus_too_long, other.failed_assigned_counts.consensus_too_long);
  AddOther(failed_assigned_counts.too_many_errors, other.failed_assigned_counts.too_many_errors);
  AddOther(failed_assigned_counts.trimmed_read_too_short, other.failed_assigned_counts.trimmed_read_too_short);
  AddOther(failed_assigned_counts.failed_hairpin_stem_trim_reads,
           other.failed_assigned_counts.failed_hairpin_stem_trim_reads);
  AddOther(failed_assigned_counts.raw_bases, other.failed_assigned_counts.raw_bases);

  // merge mid adapter counts
  AddOther(midadapter_counts.found_by_global_symmetry, other.midadapter_counts.found_by_global_symmetry);
  AddOther(midadapter_counts.found_by_local_symmetry, other.midadapter_counts.found_by_local_symmetry);
  AddOther(midadapter_counts.found_by_string_compare, other.midadapter_counts.found_by_string_compare);

  // merge strand counts
  AddOther(strand_counts.fw, other.strand_counts.fw);
  AddOther(strand_counts.rv, other.strand_counts.rv);
  AddOther(strand_counts.fw_sig, other.strand_counts.fw_sig);
  AddOther(strand_counts.rv_sig, other.strand_counts.rv_sig);

  // merge passing read counts
  AddOther(passing_counts.full_duplex, other.passing_counts.full_duplex);
  AddOther(passing_counts.partial_duplex, other.passing_counts.partial_duplex);
  AddOther(passing_counts.longer_r2, other.passing_counts.longer_r2);
  AddOther(passing_counts.longer_r2_full_duplex, other.passing_counts.longer_r2_full_duplex);
  AddOther(passing_counts.both_umi, other.passing_counts.both_umi);
  AddOther(passing_counts.only_5p_umi, other.passing_counts.only_5p_umi);
  AddOther(passing_counts.only_3p_umi, other.passing_counts.only_3p_umi);
  AddOther(passing_counts.no_endadapter, other.passing_counts.no_endadapter);

  // merge passing base counts
  AddOther(passing_counts.concordant_bases, other.passing_counts.concordant_bases);
  AddOther(passing_counts.discordant_bases, other.passing_counts.discordant_bases);
  AddOther(passing_counts.simplex_bases, other.passing_counts.simplex_bases);
  AddOther(passing_counts.raw_bases, other.passing_counts.raw_bases);

  // merge read length distributions
  AddOther(full_duplex_length_distr, other.full_duplex_length_distr);
  AddOther(partial_duplex_length_distr, other.partial_duplex_length_distr);
  AddOther(endadapter_position_distr, other.endadapter_position_distr);
  AddOther(passing_length_distr, other.passing_length_distr);
  unassigned_length_distr += other.unassigned_length_distr;
  no_hairpin_length_distr += other.no_hairpin_length_distr;
  total_length_distr += other.total_length_distr;
}

// For early stopping, returns the count of the sample with the smallest count of concordant duplex bases
size_t DuplexMetrics::MinConcordDupBases() {
  PerSIDCount total(metrics_constraints::max_sid_id_index + 1, 0);
  // Sum all the concordant_duplex_bases_total per sample
  const std::function func = [&total](const DuplexMetrics& item) {
    AddOther(total, item.passing_counts.concordant_bases);
  };
  instance.ForEach(func);
  // Return the count of the sample with the minimum concordant duplex bases
  return *std::ranges::min_element(total.begin(), total.end());
}

DuplexMetrics DuplexMetrics::SumTotal() {
  DuplexMetrics total;
  const std::function func = [&total](const DuplexMetrics& item) { total.Add(item); };
  instance.ForEach(func);
  return total;
}

void DuplexMetrics::WriteMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool,
                                 const AdapterType adapter_type) const {
  fs::create_directory(params.out_dir / kMetricsDirectory);
  WriteRunMetrics(params, sid_pool, adapter_type);
  WriteSampleMetrics(params, sid_pool, adapter_type);
  WriteReadLengthDistributions(params);
  WriteSampleAssignmentMetrics(params, sid_pool);
}

static f32 Percentage(const u64 count, const u64 total) {
  f32 percentage = 0.0f;
  if (total > 0) {
    percentage = (static_cast<f32>(count) / static_cast<f32>(total)) * 100.0f;
  }
  return percentage;
}

/// Produce a 3-column run metric row where both count and percentage are NA (metric not applicable).
static std::vector<std::string> NaRunMetricRow(const std::string_view metric) {
  return {std::string(metric), std::string(kNA), std::string(kNA)};
}

/// Produce a per-SID vector filled with NA (metric not applicable for this adapter/config).
static std::vector<std::string> ToNaVector(const size_t sid_pool_size, const std::string_view metric) {
  std::vector<std::string> result;
  result.reserve(sid_pool_size + 1);
  result.emplace_back(metric);
  result.insert(result.end(), sid_pool_size, std::string(kNA));

  return result;
}

static f64 ComputeMeanReadLength(const PerSIDHistogram& passing_length_distr) {
  u64 total_count = 0;
  u64 total_read_length_sum = 0;
  for (auto& sid_hist : passing_length_distr) {
    total_count += ComputeCount(sid_hist);
    total_read_length_sum += histogram::ComputeValueSum(sid_hist);
  }
  return total_count > 0 ? static_cast<f64>(total_read_length_sum) / static_cast<f64>(total_count) : 0.0f;
}

u64 DuplexMetrics::CountAssignedSids(const SidPool& sid_pool) const {
  return static_cast<u64>(std::ranges::count_if(sid_pool, [this](const auto& sid) {
    const auto id = sid.id;
    return id < passing_counts.full_duplex.size() &&
           passing_counts.full_duplex.at(id) + passing_counts.partial_duplex.at(id) +
                   failed_assigned_counts.consensus_too_long.at(id) + failed_assigned_counts.too_many_errors.at(id) +
                   failed_assigned_counts.trimmed_read_too_short.at(id) >
               0UL;
  }));
}

void DuplexMetrics::WriteRunMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool,
                                    const AdapterType adapter_type) const {
  // compute composite metrics
  // total assigned reads = passing reads + failed assigned reads
  // failed assigned reads metrics
  const auto consensus_too_long = accumulate(failed_assigned_counts.consensus_too_long.begin(),
                                             failed_assigned_counts.consensus_too_long.end(), 0UL);
  const auto too_many_errors =
      accumulate(failed_assigned_counts.too_many_errors.begin(), failed_assigned_counts.too_many_errors.end(), 0UL);
  const auto trimmed_read_too_short = accumulate(failed_assigned_counts.trimmed_read_too_short.begin(),
                                                 failed_assigned_counts.trimmed_read_too_short.end(), 0UL);
  const auto failed_midadapter_trim_reads{accumulate(failed_assigned_counts.failed_hairpin_stem_trim_reads.begin(),
                                                     failed_assigned_counts.failed_hairpin_stem_trim_reads.end(), 0UL)};
  const auto failed_assigned =
      consensus_too_long + too_many_errors + trimmed_read_too_short + failed_midadapter_trim_reads;

  // passing reads = full duplex reads + partial duplex reads
  const auto full_duplex_reads =
      std::accumulate(passing_counts.full_duplex.begin(), passing_counts.full_duplex.end(), 0UL);
  const auto partial_duplex_reads =
      std::accumulate(passing_counts.partial_duplex.begin(), passing_counts.partial_duplex.end(), 0UL);
  const auto passing_reads = full_duplex_reads + partial_duplex_reads;

  const auto unassigned_reads =
      unassigned_counts.read_too_short + unassigned_counts.read_too_long + unassigned_counts.no_hairpin_found;

  // total assigned reads = passing reads + failed assigned reads
  const auto assigned_reads = failed_assigned + passing_reads;
  const auto total_reads = assigned_reads + unassigned_reads;

  // mid adapter finding counts
  const auto found_string(std::accumulate(midadapter_counts.found_by_string_compare.begin(),
                                          midadapter_counts.found_by_string_compare.end(), 0UL));
  const auto found_global(std::accumulate(midadapter_counts.found_by_global_symmetry.begin(),
                                          midadapter_counts.found_by_global_symmetry.end(), 0UL));
  const auto found_local(std::accumulate(midadapter_counts.found_by_local_symmetry.begin(),
                                         midadapter_counts.found_by_local_symmetry.end(), 0UL));

  // compute passing read counts
  const auto longer_r2{std::accumulate(passing_counts.longer_r2.begin(), passing_counts.longer_r2.end(), 0UL)};
  const auto longer_r2_full_duplex{
      std::accumulate(passing_counts.longer_r2_full_duplex.begin(), passing_counts.longer_r2_full_duplex.end(), 0UL)};
  const auto both_umi{std::accumulate(passing_counts.both_umi.begin(), passing_counts.both_umi.end(), 0UL)};
  const auto only_5p_umi{std::accumulate(passing_counts.only_5p_umi.begin(), passing_counts.only_5p_umi.end(), 0UL)};
  const auto only_3p_umi{std::accumulate(passing_counts.only_3p_umi.begin(), passing_counts.only_3p_umi.end(), 0UL)};
  const auto no_endadapter{accumulate(passing_counts.no_endadapter.begin(), passing_counts.no_endadapter.end(), 0UL)};

  // compute strand counts
  const auto strand_fw_reads{std::accumulate(strand_counts.fw.begin(), strand_counts.fw.end(), 0UL)};
  const auto strand_rv_reads{std::accumulate(strand_counts.rv.begin(), strand_counts.rv.end(), 0UL)};
  const auto strand_fw_sig_reads{std::accumulate(strand_counts.fw_sig.begin(), strand_counts.fw_sig.end(), 0UL)};
  const auto strand_rv_sig_reads{std::accumulate(strand_counts.rv_sig.begin(), strand_counts.rv_sig.end(), 0UL)};

  // compute base counts
  // concordant bases = concordant duplex bases + discordant bases + non-duplex bases
  const auto concordant{
      std::accumulate(passing_counts.concordant_bases.begin(), passing_counts.concordant_bases.end(), 0UL)};
  const auto discordant{
      std::accumulate(passing_counts.discordant_bases.begin(), passing_counts.discordant_bases.end(), 0UL)};
  const auto simplex{std::accumulate(passing_counts.simplex_bases.begin(), passing_counts.simplex_bases.end(), 0UL)};
  const auto total_bases = concordant + discordant + simplex;

  const auto passing_input_bases{
      std::accumulate(passing_counts.raw_bases.begin(), passing_counts.raw_bases.end(), 0UL)};
  const auto trimmed_bases = passing_input_bases - total_bases;

  const auto failed_assigned_bases_total{
      std::accumulate(failed_assigned_counts.raw_bases.begin(), failed_assigned_counts.raw_bases.end(), 0UL)};
  const auto raw_bases = passing_input_bases + failed_assigned_bases_total + unassigned_counts.raw_bases;

  const auto out_dir = params.out_dir / kMetricsDirectory;
  const auto out_file_name{out_dir / kRunMetricsFile};

  std::ofstream out_file(out_file_name);

  // Set floating point precision for the entire file stream
  out_file << std::fixed << std::setprecision(2);

  // Create a TSV writer from the file stream
  auto writer = csv::make_tsv_writer(out_file);
  io::WriteTsvMetadata(out_file, params.command_line_info);

  // Write the header row
  writer << std::vector<std::string>{"metrics_name", "count", "percentage"};

  // Write data rows using std::make_tuple
  writer << std::make_tuple(kTotalReads, total_reads, 100.00);
  writer << std::make_tuple(kAssignedReads, assigned_reads, Percentage(assigned_reads, total_reads));

  // passing read metrics
  writer << std::make_tuple(kPassingReads, passing_reads, Percentage(passing_reads, total_reads));
  writer << std::make_tuple(kFullDuplexReads, full_duplex_reads, Percentage(full_duplex_reads, total_reads));
  writer << std::make_tuple(kPartialDuplexReads, partial_duplex_reads, Percentage(partial_duplex_reads, total_reads));

  // unassigned read metrics
  writer << std::make_tuple(kUnassignedReads, unassigned_reads, Percentage(unassigned_reads, total_reads));
  writer << std::make_tuple(kNoHairpinReads, unassigned_counts.no_hairpin_found,
                            Percentage(unassigned_counts.no_hairpin_found, total_reads));
  writer << std::make_tuple(kTooLongReads, unassigned_counts.read_too_long,
                            Percentage(unassigned_counts.read_too_long, total_reads));
  writer << std::make_tuple(kTooShortReads, unassigned_counts.read_too_short,
                            Percentage(unassigned_counts.read_too_short, total_reads));

  // failed assigned read metrics
  writer << std::make_tuple(kFailedAssignedReads, failed_assigned, Percentage(failed_assigned, total_reads));
  writer << std::make_tuple(kTooManyErrorsReads, too_many_errors, Percentage(too_many_errors, total_reads));
  writer << std::make_tuple(kTooShortTrimmedReads, trimmed_read_too_short,
                            Percentage(trimmed_read_too_short, total_reads));
  writer << std::make_tuple(kTooLongConsensusReads, consensus_too_long, Percentage(consensus_too_long, total_reads));
  writer << std::make_tuple(kFailedHairpinStemTrimReads, failed_midadapter_trim_reads,
                            Percentage(failed_midadapter_trim_reads, total_reads));

  // mid adapter finding metrics
  writer << std::make_tuple(kFoundByStringCompare, found_string, Percentage(found_string, assigned_reads));
  writer << std::make_tuple(kFoundByGlobalSymmetry, found_global, Percentage(found_global, assigned_reads));
  writer << std::make_tuple(kFoundByLocalSymmetry, found_local, Percentage(found_local, assigned_reads));

  // passing read metrics do not contribute to passing read counts (can co-occur)
  writer << std::make_tuple(kLongerR2Reads, longer_r2, Percentage(longer_r2, passing_reads));
  writer << std::make_tuple(kLongerR2FullDuplexReads, longer_r2_full_duplex,
                            Percentage(longer_r2_full_duplex, passing_reads));
  // UMI metrics — only applicable for kDuplexUMI adapters
  if (adapter_type == AdapterType::kDuplexUMI) {
    writer << std::make_tuple(kBothUmiReads, both_umi, Percentage(both_umi, passing_reads));
    writer << std::make_tuple(k5pUmiReads, only_5p_umi, Percentage(only_5p_umi, passing_reads));
    writer << std::make_tuple(k3pUmiReads, only_3p_umi, Percentage(only_3p_umi, passing_reads));
  } else {
    writer << NaRunMetricRow(kBothUmiReads);
    writer << NaRunMetricRow(k5pUmiReads);
    writer << NaRunMetricRow(k3pUmiReads);
  }
  writer << std::make_tuple(kNoEndadapterReads, no_endadapter, Percentage(no_endadapter, passing_reads));

  // strand metrics — only applicable when strand detection is enabled
  if (params.strand_detector.has_value()) {
    writer << std::make_tuple(kStrandFwReads, strand_fw_reads, Percentage(strand_fw_reads, passing_reads));
    writer << std::make_tuple(kStrandRvReads, strand_rv_reads, Percentage(strand_rv_reads, passing_reads));
    writer << std::make_tuple(kStrandFwSigReads, strand_fw_sig_reads, Percentage(strand_fw_sig_reads, passing_reads));
    writer << std::make_tuple(kStrandRvSigReads, strand_rv_sig_reads, Percentage(strand_rv_sig_reads, passing_reads));
  } else {
    writer << NaRunMetricRow(kStrandFwReads);
    writer << NaRunMetricRow(kStrandRvReads);
    writer << NaRunMetricRow(kStrandFwSigReads);
    writer << NaRunMetricRow(kStrandRvSigReads);
  }

  // base counts
  writer << std::make_tuple(kTotalBases, total_bases, 100.00);
  writer << std::make_tuple(kConcordantDuplexBases, concordant, Percentage(concordant, total_bases));
  writer << std::make_tuple(kDiscordantDuplexBases, discordant, Percentage(discordant, total_bases));
  writer << std::make_tuple(kDuplexBases, concordant + discordant, Percentage(concordant + discordant, total_bases));
  writer << std::make_tuple(kNonDuplexBases, simplex, Percentage(simplex, total_bases));

  writer << std::make_tuple(kRawBases, raw_bases, 100.00);
  writer << std::make_tuple(kUnassignedBases, unassigned_counts.raw_bases,
                            Percentage(unassigned_counts.raw_bases, raw_bases));
  writer << std::make_tuple(kFailedAssignedBases, failed_assigned_bases_total,
                            Percentage(failed_assigned_bases_total, raw_bases));
  writer << std::make_tuple(kTrimmedBases, trimmed_bases, Percentage(trimmed_bases, raw_bases));

  // sids
  writer << std::make_tuple(kNumExpectedSids, sid_pool.size(), 100.00);

  const auto num_assigned_sids = CountAssignedSids(sid_pool);
  writer << std::make_tuple(kNumSids, num_assigned_sids, Percentage(num_assigned_sids, sid_pool.size()));

  auto mean_read_length = ComputeMeanReadLength(passing_length_distr);
  std::string mean_read_length_str = fmt::format("{:.{}f}", mean_read_length, histogram::kMeanPrecision);
  writer << std::make_tuple(kMeanPassingReadLength, mean_read_length_str, "100.00");
  out_file.close();
}

static std::vector<std::string> ToStringVector(const SidPool& sid_pool, const std::string_view metric,
                                               const PerSIDCount& count) {
  std::vector<std::string> result;
  result.reserve(sid_pool.size() + 1);
  result.emplace_back(metric);
  for (const auto& sid : sid_pool) {
    if (sid.id < count.size()) {
      result.emplace_back(std::to_string(count.at(sid.id)));
    } else {
      result.emplace_back("0");
    }
  }
  return result;
}

template <typename ComputeFunc>
static std::vector<std::string> ComputePerSidMetric(const SidPool& sid_pool, const std::string_view metric,
                                                    const ComputeFunc& compute_value) {
  std::vector<std::string> result;
  result.reserve(sid_pool.size() + 1);
  result.emplace_back(metric);
  for (const auto& sid : sid_pool) {
    const auto value = compute_value(sid.id);
    if constexpr (std::is_same_v<std::decay_t<decltype(value)>, std::string>) {
      result.emplace_back(value);
    } else {
      result.emplace_back(std::to_string(value));
    }
  }
  return result;
}

void DuplexMetrics::WriteSampleMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool,
                                       const AdapterType adapter_type) const {
  const auto out_file_name{params.out_dir / kMetricsDirectory / kSampleMetricsFile};

  std::ofstream out_file(out_file_name);
  // Create a TSV writer from the file stream
  auto writer = csv::make_tsv_writer(out_file);
  io::WriteTsvMetadata(out_file, params.command_line_info);

  // write the header row
  std::vector<std::string> current_row{kMetric};
  current_row.reserve(sid_pool.size() + 1);
  for (const auto& sid : sid_pool) {
    current_row.emplace_back(sid.name);
  }
  writer << current_row;

  // index sequence row
  current_row.assign(1, kIndexSequence);
  for (const auto& sid : sid_pool) {
    current_row.emplace_back(sid.sequence);
  }
  writer << current_row;

  // assigned reads
  writer << ComputePerSidMetric(sid_pool, kAssignedReads, [this](const u32 id) -> u64 {
    if (id < passing_counts.full_duplex.size()) {
      return passing_counts.full_duplex.at(id) + passing_counts.partial_duplex.at(id) +
             failed_assigned_counts.consensus_too_long.at(id) + failed_assigned_counts.too_many_errors.at(id) +
             failed_assigned_counts.trimmed_read_too_short.at(id) +
             failed_assigned_counts.failed_hairpin_stem_trim_reads.at(id);
    }
    return 0;
  });

  // passing reads and related counts
  writer << ComputePerSidMetric(sid_pool, kPassingReads, [this](const u32 id) -> u64 {
    if (id < passing_counts.full_duplex.size()) {
      return passing_counts.full_duplex.at(id) + passing_counts.partial_duplex.at(id);
    }
    return 0;
  });

  writer << ToStringVector(sid_pool, kFullDuplexReads, passing_counts.full_duplex);
  writer << ToStringVector(sid_pool, kPartialDuplexReads, passing_counts.partial_duplex);

  // failed assigned reads
  writer << ComputePerSidMetric(sid_pool, kFailedAssignedReads, [this](const u32 id) -> u64 {
    if (id < failed_assigned_counts.consensus_too_long.size()) {
      return failed_assigned_counts.consensus_too_long.at(id) + failed_assigned_counts.too_many_errors.at(id) +
             failed_assigned_counts.trimmed_read_too_short.at(id) +
             failed_assigned_counts.failed_hairpin_stem_trim_reads.at(id);
    }
    return 0;
  });

  // failed assigned read counts
  writer << ToStringVector(sid_pool, kTooManyErrorsReads, failed_assigned_counts.too_many_errors);
  writer << ToStringVector(sid_pool, kTooShortTrimmedReads, failed_assigned_counts.trimmed_read_too_short);
  writer << ToStringVector(sid_pool, kTooLongConsensusReads, failed_assigned_counts.consensus_too_long);
  writer << ToStringVector(sid_pool, kFailedHairpinStemTrimReads,
                           failed_assigned_counts.failed_hairpin_stem_trim_reads);

  // midadapter finding method counts
  writer << ToStringVector(sid_pool, kFoundByStringCompare, midadapter_counts.found_by_string_compare);
  writer << ToStringVector(sid_pool, kFoundByGlobalSymmetry, midadapter_counts.found_by_global_symmetry);
  writer << ToStringVector(sid_pool, kFoundByLocalSymmetry, midadapter_counts.found_by_local_symmetry);

  // passing reads with unique properties
  writer << ToStringVector(sid_pool, kLongerR2Reads, passing_counts.longer_r2);
  writer << ToStringVector(sid_pool, kLongerR2FullDuplexReads, passing_counts.longer_r2_full_duplex);
  // UMI metrics — only applicable for kDuplexUMI adapters
  if (adapter_type == AdapterType::kDuplexUMI) {
    writer << ToStringVector(sid_pool, kBothUmiReads, passing_counts.both_umi);
    writer << ToStringVector(sid_pool, k5pUmiReads, passing_counts.only_5p_umi);
    writer << ToStringVector(sid_pool, k3pUmiReads, passing_counts.only_3p_umi);
  } else {
    writer << ToNaVector(sid_pool.size(), kBothUmiReads);
    writer << ToNaVector(sid_pool.size(), k5pUmiReads);
    writer << ToNaVector(sid_pool.size(), k3pUmiReads);
  }
  writer << ToStringVector(sid_pool, kNoEndadapterReads, passing_counts.no_endadapter);

  // strand metrics — only applicable when strand detection is enabled
  if (params.strand_detector.has_value()) {
    writer << ToStringVector(sid_pool, kStrandFwReads, strand_counts.fw);
    writer << ToStringVector(sid_pool, kStrandRvReads, strand_counts.rv);
    writer << ToStringVector(sid_pool, kStrandFwSigReads, strand_counts.fw_sig);
    writer << ToStringVector(sid_pool, kStrandRvSigReads, strand_counts.rv_sig);
  } else {
    writer << ToNaVector(sid_pool.size(), kStrandFwReads);
    writer << ToNaVector(sid_pool.size(), kStrandRvReads);
    writer << ToNaVector(sid_pool.size(), kStrandFwSigReads);
    writer << ToNaVector(sid_pool.size(), kStrandRvSigReads);
  }

  // base counts
  writer << ComputePerSidMetric(sid_pool, kTotalBases, [this](const u32 id) -> u64 {
    if (id < passing_counts.concordant_bases.size()) {
      return passing_counts.concordant_bases.at(id) + passing_counts.discordant_bases.at(id) +
             passing_counts.simplex_bases.at(id);
    }
    return 0;
  });

  writer << ToStringVector(sid_pool, kConcordantDuplexBases, passing_counts.concordant_bases);
  writer << ToStringVector(sid_pool, kDiscordantDuplexBases, passing_counts.discordant_bases);

  // create a PerSid Count for duplex bases by summing concordant and discordant counts
  writer << ComputePerSidMetric(sid_pool, kDuplexBases, [this](const u32 id) -> u64 {
    if (id < passing_counts.concordant_bases.size() && id < passing_counts.discordant_bases.size()) {
      return passing_counts.concordant_bases.at(id) + passing_counts.discordant_bases.at(id);
    }
    return 0;
  });

  writer << ToStringVector(sid_pool, kNonDuplexBases, passing_counts.simplex_bases);

  // base accounting metrics per SID
  writer << ToStringVector(sid_pool, kFailedAssignedBases, failed_assigned_counts.raw_bases);

  writer << ComputePerSidMetric(sid_pool, kTrimmedBases, [this](const u32 id) -> u64 {
    if (id < passing_counts.raw_bases.size() && id < passing_counts.concordant_bases.size() &&
        id < passing_counts.discordant_bases.size() && id < passing_counts.simplex_bases.size()) {
      const auto total_bases = passing_counts.concordant_bases.at(id) + passing_counts.discordant_bases.at(id) +
                               passing_counts.simplex_bases.at(id);
      return passing_counts.raw_bases.at(id) - total_bases;
    }
    return 0;
  });

  writer << ComputePerSidMetric(sid_pool, kMeanPassingReadLength, [this](const u32 id) -> std::string {
    const auto mean = (id < passing_length_distr.size() && !passing_length_distr.at(id).IsEmpty())
                          ? ComputeMean(passing_length_distr.at(id))
                          : 0.0;
    return fmt::format("{:.{}f}", mean, histogram::kMeanPrecision);
  });

  out_file.close();
}

static void WriteLengthHistogram(const LengthHistogram& histo, const std::string& title, const fs::path& out_file_name,
                                 const io::CommandLineInfo& command_line_info) {
  histogram::Histograms<u64> histograms;
  histograms.emplace_back(title, histo);
  histogram::WriteHistogramsToTsv(histograms, out_file_name, "length", histogram::kMaxLastBin, {}, command_line_info);
}

static void WriteSampleLengthHistogram(const PerSIDHistogram& histo, const SidPool& sid_pool,
                                       const fs::path& out_file_name, const io::CommandLineInfo& command_line_info) {
  histogram::Histograms<u64> histograms;
  histograms.reserve(sid_pool.size());
  for (const auto& sid : sid_pool) {
    histograms.emplace_back(sid.name, histo[sid.id]);
  }
  histogram::WriteHistogramsToTsv(histograms, out_file_name, "length", histogram::kMaxLastBin, {}, command_line_info);
}

void DuplexMetrics::WriteReadLengthDistributions(const DemuxAndTrimParam& params) const {
  const auto out_dir = params.out_dir / kMetricsDirectory;
  const auto& info = params.command_line_info;
  WriteLengthHistogram(total_length_distr, "total", out_dir / kTotalReadLengthDistr, info);
  WriteLengthHistogram(unassigned_length_distr, "unassigned", out_dir / kUnassignedReadLengthDistr, info);
  WriteLengthHistogram(no_hairpin_length_distr, "no_hairpin", out_dir / kNoHairpinReadLengthDistr, info);
}

void DuplexMetrics::WriteSampleAssignmentMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool) const {
  const auto out_dir = params.out_dir / kMetricsDirectory;
  const auto& info = params.command_line_info;
  WriteSampleLengthHistogram(full_duplex_length_distr, sid_pool, out_dir / kFullDuplexReadLengthDistr, info);
  WriteSampleLengthHistogram(partial_duplex_length_distr, sid_pool, out_dir / kPartialDuplexReadLengthDistr, info);
  WriteSampleLengthHistogram(endadapter_position_distr, sid_pool, out_dir / kEndadapterPositionDistr, info);
  WriteSampleLengthHistogram(passing_length_distr, sid_pool, out_dir / kPassingReadLengthDistr, info);
}

void ReportStrandMetrics(const DuplexMetrics& global_results) {
  const auto strand_fw_reads{
      std::accumulate(global_results.strand_counts.fw.begin(), global_results.strand_counts.fw.end(), 0UL)};
  const auto strand_rv_reads{
      std::accumulate(global_results.strand_counts.rv.begin(), global_results.strand_counts.rv.end(), 0UL)};
  const auto strand_fw_sig_reads{
      std::accumulate(global_results.strand_counts.fw_sig.begin(), global_results.strand_counts.fw_sig.end(), 0UL)};
  const auto strand_rv_sig_reads{
      std::accumulate(global_results.strand_counts.rv_sig.begin(), global_results.strand_counts.rv_sig.end(), 0UL)};
  const auto total_fw_reads = strand_fw_reads + strand_fw_sig_reads;
  const auto total_rv_reads = strand_rv_reads + strand_rv_sig_reads;
  const auto full_duplex_reads = std::accumulate(global_results.passing_counts.full_duplex.begin(),
                                                 global_results.passing_counts.full_duplex.end(), 0UL);
  const auto partial_duplex_reads = std::accumulate(global_results.passing_counts.partial_duplex.begin(),
                                                    global_results.passing_counts.partial_duplex.end(), 0UL);
  const auto passing_reads = full_duplex_reads + partial_duplex_reads;
  const auto ambiguous_reads = passing_reads - total_fw_reads - total_rv_reads;
  const auto fw_percentage = static_cast<double>(total_fw_reads) / static_cast<double>(passing_reads);
  const auto rv_percentage = static_cast<double>(total_rv_reads) / static_cast<double>(passing_reads);
  const auto ambiguous_reads_percentage = static_cast<double>(ambiguous_reads) / static_cast<double>(passing_reads);
  Logging::Info(
      "Strand detection: Forward reads: {} ({:.2f}%), Reverse reads: {} ({:.2f}%), "
      "Ambiguous reads: {} ({:.2f}%)",
      total_fw_reads, fw_percentage * 100.0, total_rv_reads, rv_percentage * 100.0, ambiguous_reads,
      ambiguous_reads_percentage * 100.0);
}

}  // namespace xoos::demux
