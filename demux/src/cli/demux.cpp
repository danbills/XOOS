#include "demux.h"

#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/file-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/enum/enum-util.h>
#include <xoos/util/file-functions.h>

#include <algorithm>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "io/sequence-reader.h"
#include "task/batch.h"
#include "task/flow-context.h"

namespace xoos::demux {

using cli::AddEnumOption;
using cli::AddOutputFileOption;
using cli::AddThreadCountOption;

/**
 * Return all sequence files from a directory recursive manner
 */
static std::vector<fs::path> RecursiveFindSequenceFiles(const fs::path& input) {
  std::vector<fs::path> stack{input};

  std::vector<fs::path> files;
  while (!stack.empty()) {
    const auto current = stack.back();
    stack.pop_back();

    if (fs::is_directory(current)) {
      for (const auto& file : fs::directory_iterator(current)) {
        stack.push_back(file.path());
      }
    } else if (IsSequenceFileFormat(current)) {
      files.push_back(current);
    }
  }
  return files;
}

static std::vector<fs::path> ReadInputFileList(const fs::path& input_file_list) {
  std::ifstream ifs{input_file_list};
  std::vector<fs::path> input_files;
  std::string line;
  while (std::getline(ifs, line)) {
    if (!fs::exists(line)) {
      throw CLI::ValidationError(fmt::format("Input does not exist: '{}'", line));
    }
    input_files.emplace_back(line);
  }
  return input_files;
}

std::vector<fs::path> ExpandInputFileList(const std::vector<fs::path>& input_files) {
  std::vector<fs::path> expanded_input_files;
  for (const auto& input_file : input_files) {
    const auto files{RecursiveFindSequenceFiles(input_file)};
    expanded_input_files.insert(expanded_input_files.end(), files.cbegin(), files.cend());
  }
  return expanded_input_files;
}

/**
 * @brief Validator for --adapter-design-name.
 *
 * Side effect: enables methylation mode when the SBX-DM design is selected.
 *
 * @param[in,out] param Parameter bundle to mutate when applicable.
 * @param[in]     input User-supplied adapter design name.
 * @return Empty string on success; an error message describing the problem on failure.
 */
static std::string ValidateAdapterDesignName(DemuxAndTrimParam& param, const std::string_view input) {
  if (input.empty()) {
    return "adapter-design-name must not be empty";
  }
  if (input == "SBX-DM") {
    param.methylation = true;
  }
  return {};
}

/**
 * @brief Configure the strand detector and classifier from a precomputed strand index file.
 *
 * Loads the index, logs the per-strand false-positive rates, and constructs a strand classifier
 * using the average FPR and the (already-parsed) Phred-scaled critical value.
 *
 * @param[in,out] param    Parameter bundle holding the strand_critical value and receiving the
 *                         constructed strand_detector / strand_classifier.
 * @param[in]     filename Path to the strand index file produced by demux-strand-index.
 */
static void ConfigureStrandDetector(DemuxAndTrimParam& param, const fs::path& filename) {
  param.strand_detector.emplace(filename);
  const auto [fw_fpr, rv_fpr] = param.strand_detector->FalsePositiveRates();
  Logging::Info("Strand Detection Mode. Filter FPRs fw: {:.5f} rv: {:.5f}", fw_fpr, rv_fpr);
  const f64 average_error_rate = (fw_fpr + rv_fpr) / 2.0;
  const f64 rescaled_critical_value = std::pow(10, -param.strand_critical / 10.0);
  param.strand_classifier.emplace(*param.strand_detector, average_error_rate, rescaled_critical_value);
}

/**
 * @brief Define CLI options specific to simplex adapter processing.
 *
 * @param[in]     app   CLI11 application to register options with.
 * @param[in,out] param Shared parameter bundle that receives parsed values.
 */
static void DefineSimplexOptions(CLI::App* const app, const std::shared_ptr<DemuxAndTrimParam>& param) {
  AddEnumOption(
      app, "--read-length-mode", param->read_length_mode,
      fmt::format("Simplex read output mode. "
                  "{}: emit full and partial reads together. "
                  "{}: emit only full-length reads. "
                  "{}: emit full and partial reads into separate subdirectories (<sample>/{}/, <sample>/{}).",
                  enum_util::FormatEnumName(ReadLengthMode::kAll), enum_util::FormatEnumName(ReadLengthMode::kFullOnly),
                  enum_util::FormatEnumName(ReadLengthMode::kAllSplit), kFullReadSubdir, kPartialReadSubdir),
      param->read_length_mode)
      ->group("Simplex Options");
  AddEnumOption(app, "--discordant-sid-mode", param->discordant_sid_mode,
                fmt::format("Determines how to handle reads where 5' and 3' SIDs disagree. "
                            "{}: discard only when both sides are equally likely (tied). "
                            "{}: discard all discordant reads (both sides significantly detected). "
                            "{}: keep all discordant reads, assign to winning side. Ties are set to 5'.",
                            enum_util::FormatEnumName(DiscordantSidMode::kDiscardTied),
                            enum_util::FormatEnumName(DiscordantSidMode::kDiscardAll),
                            enum_util::FormatEnumName(DiscordantSidMode::kKeep)),
                param->discordant_sid_mode)
      ->group("Simplex Options");
  app->add_flag(
         "--suppress-simplex-qual-override", param->suppress_simplex_qual_override,
         "Suppress overriding base-call quality scores for simplex adapter reads (YS, YSU, YS-NEW, SIMPLEX-10X). "
         "By default, the scores are replaced with a fixed simplex value. "
         "When this flag is set, the original base-call quality scores are preserved.")
      ->default_val(param->suppress_simplex_qual_override)
      ->group("Simplex Options");
  app->add_option("--min-score", param->min_score, "Minimum log odds based score needed for a valid adapter match.")
      ->check(CLI::Range(s32{0}, std::numeric_limits<s32>::max()))
      ->default_val(param->min_score)
      ->group("Simplex Options");
}

/**
 * @brief Define CLI options specific to duplex adapter processing and strand detection.
 *
 * @param[in]     app   CLI11 application to register options with.
 * @param[in,out] param Shared parameter bundle that receives parsed values.
 */
static void DefineDuplexOptions(CLI::App* const app, const std::shared_ptr<DemuxAndTrimParam>& param) {
  app->add_option("--max-error-rate-percent", param->max_error_rate_percent,
                  "Filter duplex reads with a consensus region with error rate larger than this percentage.")
      ->default_val(param->max_error_rate_percent)
      ->group("Duplex Options");
  app->add_option("--stop-after-min-concordant-duplex-bases", param->min_concord_dp_bases,
                  "Minimum concordant duplex bases per sample for early stopping.")
      ->group("Duplex Options");
  app->add_option("--strand-critical-phred", param->strand_critical,
                  "In duplex strand detection mode, this is the Phred scaled FPR critical value "
                  "-10*log10(critical_value). Low values will classify faster at cost of accuracy.")
      ->default_val(param->strand_critical)
      ->group("Duplex Options");
  app->add_option_function<fs::path>(
         "--strand-detect", [&param](const fs::path& filename) { ConfigureStrandDetector(*param, filename); },
         "Enable genome strand detection with the specified file created by demux_strand_index.")
      ->group("Duplex Options");
}

void DefineOptions(CLI::App* const app, const std::shared_ptr<DemuxAndTrimParam>& param) {
  app->add_option_function<std::vector<fs::path>>(
         "-i,--input",
         [&param](const std::vector<fs::path>& inputs) {
           const auto expanded{ExpandInputFileList(inputs)};
           std::ranges::copy(expanded, std::back_inserter(param->inputs));
         },
         "Inputs containing untrimmed sequences in either FASTQ, FASTQ.gz, or RDB format. "
         "May specify a directory (reading all relevant files in directory). May be specified multiple times, or as a "
         "space separated list.")
      ->check(CLI::ExistingPath);
  app->add_option_function<fs::path>(
         "--input-files-list",
         [&param](const std::string& input_file_list) {
           const auto inputs = ExpandInputFileList(ReadInputFileList(input_file_list));
           std::ranges::copy(inputs, std::back_inserter(param->inputs));
         },
         "File listing the input file paths instead of setting them via -i/--input option.")
      ->check(CLI::ExistingFile);
  app->add_option("-o,--out-dir", param->out_dir, "Output directory")->default_val(param->out_dir);
  app->add_flag("--overwrite", param->overwrite, "Overwrites existing non-empty output directory and its contents.")
      ->default_val(param->overwrite);
  app->add_option("-b,--batch-size", param->batch_size, "Number of reads to process in each batch.")
      ->default_val(param->batch_size);
  AddThreadCountOption(app, "--threads", param->threads);
  app->add_option("-p,--adapter-design-bundle", param->adapter_design_bundle,
                  "Path to ZIP file or directory containing the adapter architecture, search strategy, and default "
                  "adapter design")
      ->default_val(param->adapter_design_bundle)
      ->check(CLI::ExistingPath);
  app->add_option("-n,--adapter-design-name", param->adapter_design_name,
                  "The adapter design name to be loaded from the bundle. Defaults to SBX-D.")
      ->check(CLI::Validator([&param](const std::string& input) { return ValidateAdapterDesignName(*param, input); },
                             "valid-adapter-design-name"))
      ->default_val(param->adapter_design_name);
  app->add_option(
         "--sample-sheet", param->sample_sheet,
         "If specified only look for the adapter SIDs in sample sheet. Otherwise adapter design bundle sheet is used.")
      ->check(CLI::ExistingFile);
  app->add_option(
         "--compression-level", param->compression_level,
         "Compression level for output FASTQ if using gzip (1-9) or zstd (1-19) compression. Higher level increases "
         "compression at the cost of speed.  If compression type is set to none, this option is ignored and "
         "output will be uncompressed.  Will choose default for both algorithms if not specified.")
      ->check(CLI::Range(kMinCompressionLevel, kMaxZstdCompressionLevel));
  AddThreadCountOption(app, "--writing-threads-per-sample", param->writing_threads_per_sample)
      ->description("Number of writing threads used per sample.")
      ->check(CLI::Range(size_t{1}, BatchTask::kMaxNumberOfSinkWorkers))
      ->default_val(param->writing_threads_per_sample);
  app->add_option("--worker-threads-per-input", param->worker_threads_per_input,
                  "Controls how many input files are processed concurrently: "
                  "concurrent inputs = threads / worker-threads-per-input (minimum 1). "
                  "Higher values reduce concurrency, which can lower memory pressure for large files.")
      ->check(CLI::Range(size_t{1}, std::numeric_limits<size_t>::max()))
      ->default_val(param->worker_threads_per_input);
  AddEnumOption(
      app, "--compression-type", param->compression_type,
      "Compression type for output FASTQs. If set to none, output will ignore the compression level parameter.",
      param->compression_type);
  app->add_option("--length-distribution-report-max", param->length_distribution_report_max,
                  "The maximum size the largest bucket in the length histogram distribution report can be.")
      ->default_val(param->length_distribution_report_max);
  app->add_option("--min-read-len", param->min_read_len, "Filter raw reads shorter than this length.")
      ->check(CLI::PositiveNumber)
      ->default_val(param->min_read_len);
  app->add_option("--min-trimmed-read-len", param->min_trimmed_read_len,
                  "Filter trimmed reads shorter than this length.")
      ->check(CLI::PositiveNumber)
      ->default_val(param->min_trimmed_read_len);
  app->add_flag("--output-failed-reads", param->output_failed_reads,
                fmt::format("Output unaltered failed reads into the '{}' subdirectory under the output directory. "
                            "Files are named with prefix '{}'.",
                            kReservedRawFailedName, kReservedRawFailedName))
      ->default_val(param->output_failed_reads);
  // Hidden option for emitting trim coordinates in output FASTQ comments for trimming evaluation
  app->add_flag("--trimming-debug", param->trimming_debug)->default_val(param->trimming_debug)->group("");

  DefineSimplexOptions(app, param);
  DefineDuplexOptions(app, param);
}

static void ValidateOptions(const std::shared_ptr<DemuxAndTrimParam>& param) {
  // Make sure that we have at least one input file.
  if (param->inputs.empty()) {
    throw CLI::ValidationError("Must provide at least one input");
  }
  // Sort inputs for deterministic output filename assignment.
  std::ranges::sort(param->inputs);
  // Guard against hardware_concurrency() returning 0 when --threads 0 is used.
  if (param->threads == 0) {
    Logging::Warn("--threads resolved to 0 (hardware_concurrency unavailable). Defaulting to 1.");
    param->threads = 1;
  }
  // Clamp worker_threads_per_input to the thread count so the default (4) works with any --threads value.
  if (param->worker_threads_per_input > param->threads) {
    Logging::Warn("--worker-threads-per-input ({}) exceeds --threads ({}). Clamping to {}.",
                  param->worker_threads_per_input, param->threads, param->threads);
    param->worker_threads_per_input = param->threads;
  }
  // Make sure that min_trimmed_read_len less than or equal to min_read_len.
  // It is impossible trimmed reads to be longer than raw reads.
  if (param->min_trimmed_read_len > param->min_read_len) {
    Logging::Warn(
        "--min-trimmed-read-len ({}) is greater than --min-read-len ({}). "
        "Trimmed reads cannot be longer than raw reads so we will set --min-trimmed-read-len to --min-read-len.",
        param->min_trimmed_read_len, param->min_read_len);
    param->min_trimmed_read_len = param->min_read_len;
  }

  // if the output directory is a file, then throw an error
  if (fs::exists(param->out_dir) && !fs::is_directory(param->out_dir)) {
    throw CLI::ValidationError(fmt::format("Output path ({}) is not a directory.", param->out_dir));
  }

  // create directory if it does not exist
  if (!fs::exists(param->out_dir)) {
    fs::create_directories(param->out_dir);
  }
  // if it exists, check if we can write to it
  if (fs::exists(param->out_dir)) {
    // throw an error if unable to write to the directory
    file::CheckFileIsWritable(param->out_dir);
  } else {
    // not created so throwing an error
    throw CLI::ValidationError(
        fmt::format("Output directory ({}) does not exist and could not be created.", param->out_dir));
  }
}

static void PreMainProcessing(const cli::ConstAppPtr app, const std::shared_ptr<DemuxAndTrimParam>& param) {
  ValidateOptions(param);
  param->command_line_info = cli::GetCommandLineInfo(app);
}

cli::PreCallback<DemuxAndTrimParam> CreatePreCallback() {
  return [](const cli::ConstAppPtr app, const std::shared_ptr<DemuxAndTrimParam>& param) {
    PreMainProcessing(app, param);
  };
}

}  // namespace xoos::demux
