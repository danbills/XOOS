#include "generate-pon.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <queue>
#include <ranges>
#include <set>
#include <string>
#include <string_view>
#include <thread>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-writer.h>
#include <xoos/log/logging.h>
#include <xoos/util/string-functions.h>

#include "util/file-util.h"

namespace xoos::svc {

// PON-specific VCF INFO field names
static constexpr std::string_view kPonOcc = "PON_OCC";
static constexpr std::string_view kPonOccRatio = "PON_OCC_RATIO";
static constexpr std::string_view kPonMeanDuplexAf = "PON_MEAN_DUPLEX_AF";
static constexpr std::string_view kPonMinDuplexAf = "PON_MIN_DUPLEX_AF";
static constexpr std::string_view kPonMaxDuplexAf = "PON_MAX_DUPLEX_AF";

// Expected column names in the feature TSV files
static constexpr std::string_view kColChrom = "chrom";
static constexpr std::string_view kColPos = "pos";
static constexpr std::string_view kColRef = "ref";
static constexpr std::string_view kColAlt = "alt";
static constexpr std::string_view kColDuplexAf = "duplex_af";
static constexpr std::string_view kColMapqMean = "mapq_mean";

// Indices of required columns in the feature TSV header
struct ColumnIndices {
  size_t chrom{};
  size_t pos{};
  size_t ref{};
  size_t alt{};
  size_t duplex_af{};
  size_t mapq_mean{};
};

// ── ShardedVariantMap implementation ──

void ShardedVariantMap::InsertOrUpdate(const VariantId& key, const f32 duplex_af) {
  const auto shard_idx = std::hash<VariantId>{}(key) % kNumShards;
  auto& shard = shards[shard_idx];
  const std::lock_guard lock(shard.mutex);

  const auto [it, inserted] = shard.map.try_emplace(
      key,
      VariantStats{
          .occurrence = 1U, .min_duplex_af = duplex_af, .max_duplex_af = duplex_af, .sum_duplex_af = duplex_af});
  if (!inserted) {
    auto& stats = it->second;
    ++stats.occurrence;
    stats.min_duplex_af = std::min(stats.min_duplex_af, duplex_af);
    stats.max_duplex_af = std::max(stats.max_duplex_af, duplex_af);
    stats.sum_duplex_af += duplex_af;
  }
}

size_t ShardedVariantMap::Size() const {
  size_t total = 0;
  for (const auto& shard : shards) {
    total += shard.map.size();
  }
  return total;
}

std::set<std::string, std::less<>> ShardedVariantMap::UniqueContigs() const {
  std::set<std::string, std::less<>> contigs;
  for (const auto& shard : shards) {
    for (const auto& [key, _] : shard.map) {
      contigs.insert(key.chrom);
    }
  }
  return contigs;
}

/**
 * @brief Parse the TSV header line and return column indices for required fields.
 * @param header_line Tab-separated header line from the feature file.
 * @param file_path Path to the feature file (for error messages).
 * @return Column indices for chrom, pos, ref, alt, duplex_af, mapq_mean.
 * @throws xoos::Error if any required column is missing.
 */
static ColumnIndices ParseHeader(const std::string& header_line, const fs::path& file_path) {
  const auto columns = string::Split(header_line, "\t");
  ColumnIndices indices{};

  struct ColumnMapping {
    std::string_view name;
    size_t* index;
    bool found;
  };

  std::array<ColumnMapping, 6> mappings = {{
      {kColChrom, &indices.chrom, false},
      {kColPos, &indices.pos, false},
      {kColRef, &indices.ref, false},
      {kColAlt, &indices.alt, false},
      {kColDuplexAf, &indices.duplex_af, false},
      {kColMapqMean, &indices.mapq_mean, false},
  }};

  for (size_t i = 0; i < columns.size(); ++i) {
    for (auto& [name, index, found] : mappings) {
      if (columns[i] == name) {
        *index = i;
        found = true;
        break;
      }
    }
  }

  const bool all_found = std::ranges::all_of(mappings, [](const auto& m) { return m.found; });
  if (!all_found) {
    throw error::Error(
        "Feature file {} is missing required columns. Expected: chrom, pos, ref, alt, duplex_af, "
        "mapq_mean",
        file_path);
  }
  return indices;
}

/**
 * @brief Read a single feature TSV file, filter variants, and insert into the sharded map.
 * @param file_path Path to the feature TSV file produced by compute_bam_features.
 * @param min_duplex_af Minimum duplex_af threshold (inclusive).
 * @param max_duplex_af Maximum duplex_af threshold (inclusive).
 * @param min_mapq_mean Minimum mapq_mean threshold (inclusive).
 * @param sharded_map Shared map to insert passing variants into.
 * @throws xoos::Error if the file cannot be opened or has a malformed header.
 */
static void ProcessFeatureFile(const fs::path& file_path,
                               const f32 min_duplex_af,
                               const f32 max_duplex_af,
                               const f32 min_mapq_mean,
                               ShardedVariantMap& sharded_map) {
  std::ifstream input(file_path);
  if (!input.is_open()) {
    throw error::Error("Cannot open feature file: {}", file_path);
  }

  // Find the header line (skip comment lines starting with #)
  std::string line;
  bool header_found = false;
  while (std::getline(input, line)) {
    if (line.empty() || line.starts_with(io::kTsvCommentLinePrefix)) {
      continue;
    }
    header_found = true;
    break;
  }
  if (!header_found) {
    throw error::Error("Cannot find header in feature file: {}", file_path);
  }

  string::Trim(line);
  const auto indices = ParseHeader(line, file_path);
  const auto min_columns =
      std::max({indices.chrom, indices.pos, indices.ref, indices.alt, indices.duplex_af, indices.mapq_mean}) + 1U;

  u64 total_variants = 0;
  u64 passed_variants = 0;

  while (std::getline(input, line)) {
    if (line.empty() || line.starts_with(io::kTsvCommentLinePrefix)) {
      continue;
    }
    string::Trim(line);
    const auto fields = string::Split(line, "\t");
    if (fields.size() < min_columns) {
      Logging::Warn("Skipping malformed line {} in {} (expected {} columns, got {})",
                    total_variants + 1U,
                    file_path,
                    min_columns,
                    fields.size());
      continue;
    }

    ++total_variants;

    f32 duplex_af{};
    f32 mapq_mean{};
    u64 tsv_pos{};
    try {
      duplex_af = std::stof(fields[indices.duplex_af]);
      mapq_mean = std::stof(fields[indices.mapq_mean]);
      tsv_pos = std::stoull(fields[indices.pos]);
    } catch (const std::exception& e) {
      throw error::Error("Failed to parse line {} in {}: {}", total_variants + 1U, file_path, e.what());
    }

    // Apply per-variant filters
    if (duplex_af < min_duplex_af || duplex_af > max_duplex_af || mapq_mean < min_mapq_mean) {
      continue;
    }

    ++passed_variants;

    // Feature TSV positions are 1-based (VariantInfoSerializer writes vid.pos + 1).
    // Convert to 0-based for htslib's SetPosition().
    if (tsv_pos == 0U) {
      throw error::Error("Invalid position 0 at line {} in {} (expected 1-based)", total_variants + 1U, file_path);
    }
    const VariantId key(fields[indices.chrom], tsv_pos - 1U, fields[indices.ref], fields[indices.alt]);

    sharded_map.InsertOrUpdate(key, duplex_af);
  }

  Logging::Info("Processed {}: {} total variants, {} passed filters",
                file_path.filename().string(),
                total_variants,
                passed_variants);
}

/**
 * @brief Build a VCF header with contig lines and PON INFO field definitions.
 *
 * When `command_line` is provided, a `##RocheCommandLine` provenance line is
 * prepended to the header before the contig and INFO entries.
 *
 * @param contigs Sorted set of contig names observed in the variant data.
 * @param command_line Optional CLI metadata (program name, version, arguments);
 *                     when set, emits a ##RocheCommandLine header line.
 * @return Shared pointer to the constructed VCF header.
 */
static io::VcfHeaderPtr BuildPonHeader(const std::set<std::string, std::less<>>& contigs,
                                       const std::optional<io::CommandLineInfo>& command_line) {
  auto header = io::VcfHeader::Create();

  if (command_line) {
    header->AddCustomMetaDataLine(io::GetCommandInfo(*command_line));
  }

  for (const auto& contig : contigs) {
    header->AddCustomMetaDataLine(fmt::format("##contig=<ID={}>", contig));
  }

  const auto pon_occ = std::string(kPonOcc);
  const auto pon_occ_ratio = std::string(kPonOccRatio);
  const auto pon_mean_af = std::string(kPonMeanDuplexAf);
  const auto pon_min_af = std::string(kPonMinDuplexAf);
  const auto pon_max_af = std::string(kPonMaxDuplexAf);

  using enum io::FieldType;
  header->AddInfoLine({pon_occ, "Number of samples containing this variant", io::kNumberOne, kInteger});
  header->AddInfoLine({pon_occ_ratio, "Fraction of samples containing this variant", io::kNumberOne, kFloat});
  header->AddInfoLine({pon_mean_af, "Mean duplex_af across samples", io::kNumberOne, kFloat});
  header->AddInfoLine({pon_min_af, "Minimum duplex_af across samples", io::kNumberOne, kFloat});
  header->AddInfoLine({pon_max_af, "Maximum duplex_af across samples", io::kNumberOne, kFloat});

  header->Sync();
  return header;
}

/**
 * @brief Write a single variant record to the PON VCF.
 * @param writer VCF writer to emit the record.
 * @param key Variant identifier (chrom, pos, ref, alt).
 * @param stats Aggregated statistics for this variant across samples.
 * @param total_samples Total number of input samples (for occurrence ratio).
 */
static void WritePonRecord(io::VcfWriter& writer,
                           const VariantId& key,
                           const VariantStats& stats,
                           const u32 total_samples) {
  const auto occ_ratio = static_cast<f32>(stats.occurrence) / static_cast<f32>(total_samples);
  const auto mean_duplex_af = static_cast<f32>(stats.sum_duplex_af / stats.occurrence);

  const auto pon_occ = std::string(kPonOcc);
  const auto pon_occ_ratio = std::string(kPonOccRatio);
  const auto pon_mean_af = std::string(kPonMeanDuplexAf);
  const auto pon_min_af = std::string(kPonMinDuplexAf);
  const auto pon_max_af = std::string(kPonMaxDuplexAf);

  const io::VcfIntegerFields integer_info = {{pon_occ, {static_cast<int>(stats.occurrence)}}};
  const io::VcfFloatFields float_info = {{pon_occ_ratio, {occ_ratio}},
                                         {pon_mean_af, {mean_duplex_af}},
                                         {pon_min_af, {stats.min_duplex_af}},
                                         {pon_max_af, {stats.max_duplex_af}}};
  const io::VcfStringFields string_info = {};

  const io::TypedVcfFields info_fields{
      integer_info,
      float_info,
      string_info,
      std::vector<std::string>{pon_occ, pon_occ_ratio, pon_mean_af, pon_min_af, pon_max_af}};
  const io::TypedVcfFields empty_format{{}, {}, {}};

  const auto record = writer.CreateRecord(
      key.chrom, static_cast<int>(key.pos), ".", {key.ref, key.alt}, std::nullopt, "PASS", info_fields, empty_format);
  writer.WriteRecord(record);
}

void GeneratePon(const GeneratePonParam& param) {
  const auto& feature_files = param.feature_files;
  const auto requested_threads =
      param.threads == 0U ? std::max(1U, std::thread::hardware_concurrency()) : param.threads;
  const auto threads = std::min(requested_threads, feature_files.size());

  Logging::Info("Generating PON from {} feature files using {} thread(s)", feature_files.size(), threads);
  Logging::Info(
      "Filters: duplex_af in [{}, {}], mapq_mean >= {}", param.min_duplex_af, param.max_duplex_af, param.min_mapq_mean);

  // All threads insert into a single sharded map, keeping memory close to single-threaded usage
  ShardedVariantMap sharded_map;

  if (threads <= 1U) {
    for (const auto& file_path : feature_files) {
      ProcessFeatureFile(file_path, param.min_duplex_af, param.max_duplex_af, param.min_mapq_mean, sharded_map);
    }
  } else {
    tf::Executor executor(threads);
    tf::Taskflow flow;

    for (const auto& file_path : feature_files) {
      flow.emplace([&param, &sharded_map, &fp = file_path] {
        ProcessFeatureFile(fp, param.min_duplex_af, param.max_duplex_af, param.min_mapq_mean, sharded_map);
      });
    }

    executor.run(flow).get();
  }

  const auto total_variants = sharded_map.Size();
  Logging::Info("Total unique variants after filtering: {}", total_variants);

  // Write PON VCF using k-way merge across shards (no second map needed)
  CreateParentDirectoryIfNotExists(param.output_file);
  const auto header = BuildPonHeader(sharded_map.UniqueContigs(), param.command_line);
  io::VcfWriter writer(param.output_file, header);
  writer.WriteHeader();

  // K-way merge: each shard's std::map is sorted by VariantId. Use a min-heap
  // of iterators to emit variants in globally sorted order.
  using ShardIter = VariantMap::const_iterator;

  struct HeapEntry {
    const VariantId* key;
    const VariantStats* stats;
    size_t shard_idx;

    bool operator>(const HeapEntry& other) const {
      return other.key->operator<(*key);
    }
  };

  std::priority_queue<HeapEntry, vec<HeapEntry>, std::greater<>> heap;
  vec<ShardIter> iterators(kNumShards);
  for (size_t i = 0; i < kNumShards; ++i) {
    iterators[i] = sharded_map.shards[i].map.cbegin();
    if (iterators[i] != sharded_map.shards[i].map.cend()) {
      heap.emplace(HeapEntry{&iterators[i]->first, &iterators[i]->second, i});
    }
  }

  const auto total_samples = static_cast<u32>(feature_files.size());
  size_t written = 0;
  while (!heap.empty()) {
    const auto [key, stats, shard_idx] = heap.top();
    heap.pop();
    WritePonRecord(writer, *key, *stats, total_samples);
    ++written;

    ++iterators[shard_idx];
    if (iterators[shard_idx] != sharded_map.shards[shard_idx].map.cend()) {
      heap.emplace(HeapEntry{&iterators[shard_idx]->first, &iterators[shard_idx]->second, shard_idx});
    }
  }

  writer.Flush();
  Logging::Info("Wrote PON VCF with {} variants to {}", written, param.output_file.string());
}

}  // namespace xoos::svc
