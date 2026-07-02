#include "task/flow-manager.h"

#include <xoos/error/error.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>

#include <algorithm>
#include <cmath>
#include <fstream>

#include "csv.hpp"
#include "io/sample-sheet/sample-sheet.h"
#include "lut-bundle/lut-bundle.h"
#include "metrics/duplex-metrics.h"
#include "metrics/simplex-metrics.h"
#include "task/flow-context.h"
#include "utility/stop-watch.h"

namespace xoos::demux {

// NOTE: A previous ExecutorShadow struct replicated tf::Executor's internal
// memory layout to access private thread handles for CPU affinity. It was
// removed because it is ABI-fragile — any Taskflow update that changes the
// internal layout would cause this cast to read garbage memory. If thread
// affinity is needed, it must be set from within the worker thread using
// platform APIs on the current thread. tf::Executor::this_worker_id() can be
// used only to identify the calling worker and map it to a desired CPU set;
// otherwise, request an upstream API addition.

FlowManager::FlowManager(const DemuxAndTrimParam& param, const AdapterDesign& design)
    : adapter_type(design.type),
      is_duplex([this]() {
        switch (adapter_type) {
          using enum AdapterType;
          case kDuplex:
          case kDuplexUMI:
          case kDuplexStem:
            return true;
          default:
            return false;
        }
      }()),
      param(param),
      sid_pool(LoadSidPool(param, design, _lut_bundle_ysu, _lut_bundle_ys, _lut_bundle_simplex, _lut_bundle_simplex_10x,
                           _lut_bundle_duplex, _lut_bundle_duplex_umi, _lut_bundle_duplex_stem)) {
  if (is_duplex && param.read_length_mode == ReadLengthMode::kAllSplit) {
    throw error::Error("--read-length-mode all-split is only supported for simplex adapter types");
  }

  // initialize Metrics configuration
  metrics_constraints::max_sid_id_index = sid_pool.size() - 1;
  metrics_constraints::max_logged_read_length = param.length_distribution_report_max;

  InitDemuxObject(design);
  _start_time = std::chrono::high_resolution_clock::now();

  // Chunk input files so that each concurrently-processed file gets worker_threads_per_input threads.
  const size_t chunk_size = std::max(1UL, param.threads / param.worker_threads_per_input);
  Logging::Info("Demultiplexing with {} threads, {} worker threads per input, {} writer threads per sample.",
                param.threads, param.worker_threads_per_input, param.writing_threads_per_sample);
  _executor = std::make_unique<tf::Executor>(param.threads);

  for (size_t i = 0; i < param.inputs.size(); i += chunk_size) {
    for (size_t j = i; j < param.inputs.size() && j < (i + chunk_size); ++j) {
      // Creating a new flow context will add a graph to the executor.
      _flow_contexts.emplace_back(std::make_unique<FlowContext>(*this, param.inputs[j], j));
      _total_file_size += std::filesystem::file_size(param.inputs[j]);
      ++_nr_files;
    }
    // While we're processing the input, TF will keep generating new tasks for every output file. Eventually,
    // those tasks will be done and the TF will shut down. We need to wait for that to happen.
    _executor->wait_for_all();
    for (const auto& flow_context : _flow_contexts) {
      _nr_sequences += flow_context->NumSequences();
      _nr_bases += flow_context->NumBases();
      _nr_bytes += flow_context->NumBytes();
    }
    // If we arrive here, all tasks are done. We can now safely remove the graphs from the executor.
    _flow_contexts.clear();
    if (Task::HasException()) {
      std::rethrow_exception(Task::GetException());
    }
    if (param.min_concord_dp_bases.has_value()) {
      auto min_concord_dp_bases = DuplexMetrics::MinConcordDupBases();

      // Check if we need to early stop
      if (param.min_concord_dp_bases.value() < min_concord_dp_bases) {
        Logging::Info("Early stopping due to minimum concordant duplex bases reached at {} bases",
                      min_concord_dp_bases);
        break;
      }
    }
  }
  _executor.reset();
  RenameOutputFiles();
  _end_time = std::chrono::high_resolution_clock::now();
}

void FlowManager::InitDemuxObject(const AdapterDesign& design) {
  switch (design.type) {
    using enum AdapterType;
    case kYsu: {
      auto enable_partial = param.read_length_mode != ReadLengthMode::kFullOnly;
      _demux_and_trim_ysu = std::make_unique<DemuxAndTrimYsu>(enable_partial, *_lut_bundle_ysu);
      break;
    }
    case kYs: {
      auto enable_partial = param.read_length_mode != ReadLengthMode::kFullOnly;
      _demux_and_trim_ys = std::make_unique<DemuxAndTrimYs>(enable_partial, *_lut_bundle_ys);
      break;
    }
    case kSimplex:
      _demux_and_trim_simplex = std::make_unique<DemuxAndTrimSimplex>(param, *_lut_bundle_simplex);
      break;
    case kSimplex10x: {
      _demux_and_trim_simplex_10x = std::make_unique<DemuxAndTrimSimplex10x>(param, *_lut_bundle_simplex_10x);
      const auto& b = *_lut_bundle_simplex_10x;
      _demux_and_trim_simplex_10x->InitFastPath(b.fixed_1_matcher.GetFrontSeq(), b.sid_2_matcher.Pool(),
                                                b.fixed_3_matcher.GetFrontSeq());
      break;
    }
    case kDuplex:
      _hairpin_finder = std::make_unique<HairpinFinder>(*_lut_bundle_duplex);
      _demux_and_trim_duplex = std::make_unique<DemuxAndTrimDuplex>(*_lut_bundle_duplex);
      break;
    case kDuplexUMI:
      _hairpin_finder = std::make_unique<HairpinFinder>(*_lut_bundle_duplex_umi);
      _demux_and_trim_duplex_umi = std::make_unique<DemuxAndTrimDuplexUmi>(*_lut_bundle_duplex_umi);
      break;
    case kDuplexStem:
      _hairpin_finder = std::make_unique<HairpinFinder>(*_lut_bundle_duplex_stem);
      _demux_and_trim_duplex_stem = std::make_unique<DemuxAndTrimDuplexStem>(*_lut_bundle_duplex_stem);
      break;
    default:
      throw error::Error("Unsupported adapter type: {}", static_cast<s32>(design.type));
  }
}

void FlowManager::RegisterOutputPrefix(const std::string& output_prefix, const std::string_view input_filepath) const {
  std::scoped_lock lock(_provenance_mutex);
  _output_prefix_to_input[output_prefix] = input_filepath;
}

void FlowManager::RegisterOutputFile(const fs::path& output_path, const std::string& input_name_prefix) const {
  std::scoped_lock lock(_provenance_mutex);
  _output_files[output_path.parent_path()].emplace_back(output_path, input_name_prefix);
}

namespace {

// Returns the FASTQ extension including compression suffix (e.g., ".fastq.gz", ".fastq.zst", ".fastq").
std::string GetFastqExtension(const fs::path& path) {
  const auto filename = path.filename().string();
  for (const auto* ext : {".fastq.zst", ".fastq.gz", ".fastq"}) {
    if (filename.ends_with(ext)) {
      return ext;
    }
  }
  return {};
}

/**
 * @brief Derives the sample name from the output directory structure.
 *
 * If the files are written under a mode-specific subdirectory such as "<sample>/full" or
 * "<sample>/partial", the parent directory name is the sample name. Otherwise, the current
 * directory name is used directly.
 *
 * @param[in] files      Sorted list of (output_path, input_name_prefix) pairs.
 * @param[in] sample_dir Directory containing the output files for this sample.
 * @return The derived sample name.
 */
std::string DeriveSampleName(const std::vector<OutputFileEntry>& files, const fs::path& sample_dir) {
  static_cast<void>(files);

  const auto dirname = sample_dir.filename().string();
  if ((dirname == kFullReadSubdir || dirname == kPartialReadSubdir) && sample_dir.has_parent_path()) {
    // current directory name is a mode-specific subdirectory, use parent directory name as sample name
    return sample_dir.parent_path().filename().string();
  }

  // return folder name as normal
  return sample_dir.filename().string();
}

}  // namespace

void FlowManager::RenameOutputFiles() {
  constexpr size_t kDefaultPartIndexWidth = 6;

  // Upper bound on the number of output files per sample: each input file can produce up to
  // writing_threads_per_sample files.
  const auto max_file_count = param.inputs.size() * std::max(param.writing_threads_per_sample, 1UL);

  // Width is at least kDefaultPartIndexWidth, increased if needed.
  const auto part_index_width =
      std::max(kDefaultPartIndexWidth, static_cast<size_t>(std::log10(std::max(max_file_count, 1UL))) + 1);

  // Rename output files in each sample directory to the Everest naming convention
  // (<sample_name>-<part_index>.<extension>) and write a provenance TSV mapping each
  // renamed file back to its source input file.
  for (auto& [sample_dir, files] : _output_files) {
    const auto sample_name = DeriveSampleName(files, sample_dir);

    // Sort lexicographically to group files from the same input file together.
    std::ranges::sort(files);

    // Build provenance records and rename each file.
    // (new_filename, input_filepath)
    std::vector<std::pair<std::string, std::string>> provenance;

    // Assign contiguous 1-based part indices to each output file and rename it to the Everest format.
    for (size_t part_index = 1; part_index <= files.size(); ++part_index) {
      const auto& [old_path, input_name_prefix] = files[part_index - 1];
      const auto extension = GetFastqExtension(old_path);
      const auto new_filename = fmt::format("{}-{:0{}}{}", sample_name, part_index, part_index_width, extension);
      const auto new_path = sample_dir / new_filename;

      // Look up the input file path from the registered prefix.
      const auto it = _output_prefix_to_input.find(input_name_prefix);
      if (it == _output_prefix_to_input.end()) {
        throw error::Error("No input file mapping found for output prefix '{}'", input_name_prefix);
      }
      provenance.emplace_back(new_filename, it->second);

      if (fs::exists(new_path)) {
        if (param.overwrite) {
          fs::remove(new_path);
        } else {
          throw error::Error(
              "Cannot rename '{}' to '{}': target already exists. "
              "Specify --overwrite to allow overwriting existing files.",
              old_path.string(), new_path.string());
        }
      }
      fs::rename(old_path, new_path);
    }

    // Write provenance TSV in the sample directory.
    if (!provenance.empty()) {
      WriteProvenanceTsv(sample_dir, provenance);
    }
  }
}

void FlowManager::WriteProvenanceTsv(const fs::path& sample_dir,
                                     const std::vector<std::pair<std::string, std::string>>& provenance) const {
  const auto tsv_path = sample_dir / kProvenanceFilename;
  if (fs::exists(tsv_path) && !param.overwrite) {
    throw error::Error(
        "Cannot write '{}': file already exists. "
        "Specify --overwrite to allow overwriting existing files.",
        tsv_path.string());
  }
  std::ofstream tsv(tsv_path);
  io::WriteTsvMetadata(tsv, param.command_line_info);
  auto writer = csv::make_tsv_writer(tsv);
  writer << std::vector<std::string>{"output_filename", "input_filepath"};
  for (const auto& [output, input] : provenance) {
    writer << std::vector{output, input};
  }
}

// Allow SID pool to be specified in the constructor so it can exist a const object.
SidPool FlowManager::LoadSidPool(const DemuxAndTrimParam& param, const AdapterDesign& design,
                                 std::unique_ptr<LutBundleYsu>& lut_bundle_ysu,
                                 std::unique_ptr<LutBundleYs>& lut_bundle_ys,
                                 std::unique_ptr<LutBundleSimplex>& lut_bundle_simplex,
                                 std::unique_ptr<LutBundleSimplex10x>& lut_bundle_simplex_10x,
                                 std::unique_ptr<LutBundleDuplex>& lut_bundle_duplex,
                                 std::unique_ptr<LutBundleDuplexUmi>& lut_bundle_duplex_umi,
                                 std::unique_ptr<LutBundleDuplexStem>& lut_bundle_duplex_stem) {
  const auto load_bundle_sw = StopWatch{};
  auto sid_pool_list = param.sample_sheet ? std::make_optional(LoadSampleSheet(*param.sample_sheet)) : std::nullopt;
  SidPool sid_pool;
  switch (design.type) {
    using enum AdapterType;
    case kYsu:
      lut_bundle_ysu = std::make_unique<LutBundleYsu>(
          LoadLutBundle<LutBundleYsu>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = lut_bundle_ysu->Sid5pPool();
      break;
    case kYs:
      lut_bundle_ys = std::make_unique<LutBundleYs>(
          LoadLutBundle<LutBundleYs>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = lut_bundle_ys->Sid5pPool();
      break;
    case kSimplex:
      lut_bundle_simplex = std::make_unique<LutBundleSimplex>(
          LoadLutBundle<LutBundleSimplex>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = lut_bundle_simplex->sid_2_matcher.Pool();
      break;
    case kSimplex10x:
      lut_bundle_simplex_10x = std::make_unique<LutBundleSimplex10x>(
          LoadLutBundle<LutBundleSimplex10x>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = lut_bundle_simplex_10x->sid_2_matcher.Pool();
      break;
    case kDuplex:
      lut_bundle_duplex = std::make_unique<LutBundleDuplex>(
          LoadLutBundle<LutBundleDuplex>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = lut_bundle_duplex->Sid5pPool();
      break;
    case kDuplexStem:
      lut_bundle_duplex_stem = std::make_unique<LutBundleDuplexStem>(
          LoadLutBundle<LutBundleDuplexStem>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = lut_bundle_duplex_stem->Sid5pPool();
      break;
    case kDuplexUMI:
      lut_bundle_duplex_umi = std::make_unique<LutBundleDuplexUmi>(
          LoadLutBundle<LutBundleDuplexUmi>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = lut_bundle_duplex_umi->Sid5pPool();
      break;
    default:
      throw error::Error("Unsupported adapter type: {}", static_cast<s32>(design.type));
  }
  const auto load_bundle_seconds = load_bundle_sw.ElapsedTime<std::chrono::milliseconds>();
  Logging::Info("Loaded adapter design bundle in {} ms", load_bundle_seconds);

  if (sid_pool.empty()) {
    throw error::Error("SID pool is empty after loading adapter design bundle; expected at least one SID.");
  }
  // copy elision should be happening here now
  return sid_pool;
}

FlowManager::~FlowManager() {
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(_end_time - _start_time).count();
  if (is_duplex) {
    // Aggregate metrics for duplex
    const auto global_results = DuplexMetrics::SumTotal();
    global_results.WriteMetrics(param, sid_pool, adapter_type);
    // if we have strand detection turn on log info
    if (param.strand_detector.has_value()) {
      ReportStrandMetrics(global_results);
    }
  } else {
    const auto simplex_metrics_total = SimplexMetrics::SumTotal();
    simplex_metrics_total.WriteMetrics(param, sid_pool);
  }
  // Calculate Gbases per minute.
  const auto gbases = static_cast<double>(_nr_bases) * 1e-9;
  const auto gbases_per_minute = 60000.0 * gbases / static_cast<double>(duration);
  const auto mega_sequences = 1e-6 * static_cast<double>(_nr_sequences);
  const auto mega_sequences_second = 1000.0 * mega_sequences / static_cast<double>(duration);

  Logging::Info(
      "\nThroughput: processed {} sequences, {} bases, {} bytes in {} ms (using {} workers)\nNet throughput: {:.5} "
      "Gbases/minute, {:.5} MSequences/s",
      _nr_sequences, _nr_bases, _nr_bytes, duration, param.threads, gbases_per_minute, mega_sequences_second);
}

const HairpinFinder& FlowManager::HairpinFinderObj() const {
  if (!_hairpin_finder) {
    throw error::Error("FlowManager::HairpinFinderObj() called for non-duplex adapter type. Saw adapter type: {}",
                       static_cast<u32>(adapter_type));
  }
  return *_hairpin_finder;
}

const DemuxAndTrimDuplex& FlowManager::DemuxObjectDuplex() const {
  switch (adapter_type) {
    using enum AdapterType;
    case kDuplex:
      return *_demux_and_trim_duplex;
    case kDuplexStem:
      return *_demux_and_trim_duplex_stem;
    case kDuplexUMI:
      return *_demux_and_trim_duplex_umi;
    default:
      throw error::Error("FlowManager::DemuxObjectDuplex() called for non-duplex adapter type. Saw adapter type: {}",
                         static_cast<u32>(adapter_type));
  }
}

const DemuxAndTrimYsu& FlowManager::DemuxObjectYsu() const {
  if (!_demux_and_trim_ysu) {
    throw error::Error("FlowManager::DemuxObjectYsu() called but YSU demux object is not initialized");
  }
  return *_demux_and_trim_ysu;
}

const DemuxAndTrimYs& FlowManager::DemuxObjectYs() const {
  if (!_demux_and_trim_ys) {
    throw error::Error("FlowManager::DemuxObjectYs() called but YS demux object is not initialized");
  }
  return *_demux_and_trim_ys;
}

const DemuxAndTrimSimplex& FlowManager::DemuxObjectSimplex() const {
  if (!_demux_and_trim_simplex) {
    throw error::Error("FlowManager::DemuxObjectSimplex() called but simplex demux object is not initialized");
  }
  return *_demux_and_trim_simplex;
}

const DemuxAndTrimSimplex10x& FlowManager::DemuxObjectSimplex10x() const {
  if (!_demux_and_trim_simplex_10x) {
    throw error::Error("FlowManager::DemuxObjectSimplex10x() called but simplex-10x demux object is not initialized");
  }
  return *_demux_and_trim_simplex_10x;
}

}  // namespace xoos::demux
