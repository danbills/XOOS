#include "filter-variants.h"

#include <filesystem>
#include <map>
#include <memory>
#include <ranges>
#include <thread>
#include <utility>

#include <htslib/bgzf.h>
#include <htslib/vcf.h>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/io/vcf/vcf-writer.h>
#include <xoos/log/logging.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/sex_predict/vcf-sex-predictor.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/container-functions.h>
#include <xoos/util/parse-int.h>
#include <xoos/util/string-functions.h>
#include <xoos/util/tmp-dir.h>

#include "compute-bam-features/compute-bam-features.h"
#include "compute-bam-features/progress-meter.h"
#include "compute-vcf-features/compute-vcf-features.h"
#include "compute-vcf-features/vcf-header-util.h"
#include "core/filtering.h"
#include "core/model-metadata.h"
#include "core/vcf-fields.h"
#include "shap-value-tsv.h"
#include "util/bgzf-utils.h"
#include "util/file-util.h"
#include "util/log-util.h"
#include "util/parallel-compute-utils.h"

namespace xoos::svc {

using BedRegion = io::BedRegion;

/**
 * Loads a block list from file
 * @param block_list The blocklist in chr_pos_ref_alt format
 * @return The parsed left-padded blocklist
 */
static StrUnorderedSet LoadBlockList(const fs::path& block_list) {
  StrUnorderedSet result;
  std::ifstream block_list_stream(block_list);
  std::string line;
  while (getline(block_list_stream, line)) {
    // Need to left pad the blocklist so it matches the behavior of GetVariantCorrelationKey
    vec<std::string> split_line;
    std::string seg;
    std::stringstream input_line(line);
    while (std::getline(input_line, seg, '_')) {
      split_line.push_back(seg);
    }
    if (split_line.size() == 4) {
      std::string key =
          GetVariantCorrelationKey(split_line[0], util::ParseU64(split_line[1]) - 1, split_line[2], split_line[3]);
      result.insert(key);
    }
  }
  return result;
}

/**
 * @brief Helper function to extract sex chromosome name and PARs from a BED file.
 * @param bed_path Path of BED file containing PARs
 * @param default_chrom Default chromosome name if not specified in BED file
 * @return Chromosome name and vector of PARs.
 */
static std::pair<std::string, vec<Interval>> ExtractSexChromPAR(const std::optional<fs::path>& bed_path,
                                                                const std::string& default_chrom) {
  std::string name;
  vec<Interval> par;
  if (bed_path.has_value()) {
    auto chr_to_intervals = GetChromIntervalMap(bed_path).value();
    if (chr_to_intervals.empty()) {
      WarnAsErrorIfSet("No intervals found in {}", bed_path.value());
    } else if (chr_to_intervals.size() == 1) {
      name = chr_to_intervals.begin()->first;
      par = chr_to_intervals.begin()->second;
    } else {
      vec<std::string> ref_names;
      ref_names.reserve(chr_to_intervals.size());
      for (const auto& [chrom, intervals] : chr_to_intervals) {
        ref_names.emplace_back(chrom);
      }
      WarnAsErrorIfSet("PARs BED file has multiple reference names: {}", string::Join(ref_names, ","));
    }
  }
  if (par.empty()) {
    WarnAsErrorIfSet("No PARs available for chromosome {}", default_chrom);
  }
  if (name.empty()) {
    name = default_chrom;
    WarnAsErrorIfSet("Using default name for chromosome {}: {}", default_chrom, name);
  }
  return std::make_pair(name, par);
}

/**
 * @brief Extract sex chromosome information and predicted sex using pre-computed chromosomal median DP values.
 * @param param CLI parameters for `filter_variants`
 * @param chr_med_dp Pre-computed chromosome median depth values from the VCF pre-scan
 * @return Tuple: predicted sex, chrX name, chrY name, chrX PAR intervals, chrY PAR intervals, chromosomal median DPs
 */
static std::tuple<sex_predict::Sex, std::string, std::string, vec<Interval>, vec<Interval>, ChromMedianDepth>
GetSexInfo(const FilterVariantsParam& param, const ChromMedianDepth& chr_med_dp) {
  using enum sex_predict::Sex;
  // Extract the name and PARs for chromosomes X and Y
  auto [chr_x_name, chr_x_par] = ExtractSexChromPAR(param.par_bed_x, kDefaultChrXName);
  auto [chr_y_name, chr_y_par] = ExtractSexChromPAR(param.par_bed_y, kDefaultChrYName);
  auto sex = kUnknown;

  // use the "normal" sample's DP values for sex prediction
  if (chr_med_dp.normal.empty()) {
    error::Error("No median DP values extracted");
  }
  if (!chr_x_par.empty() && !chr_x_name.empty()) {
    sex = sex_predict::PredictSex(chr_med_dp.normal, param.sd_chr_name, chr_x_name);
    Logging::Info("Sex: {}", GetDescriptionForSex(sex));
    if (sex == kUnknown) {
      WarnAsErrorIfSet("Cannot determine the biological sex of input sample");
    }
  }
  return std::make_tuple(sex, chr_x_name, chr_y_name, chr_x_par, chr_y_par, chr_med_dp);
}

void FilterVariantsClass::SetGlobalContext() {
  using enum Workflow;
  // The GlobalContext struct contains various parameters required for workflow-specific filtering, and they do not
  // change between regions. The exact parameters may differ for each workflow, but certain parameters are only
  // applicable to specific workflows.

  _global_ctx.model_config = _model_config;

  // BAM feature extraction parameters
  _global_ctx.bam_feat_params =
      ComputeBamFeaturesParams{.feature_cols = _model_config.feature_cols,
                               .min_bq = _param.min_bq,
                               .min_mapq = _param.min_mapq,
                               .min_allowed_distance_from_end = _param.min_allowed_distance_from_end,
                               .min_family_size = _min_family_size,
                               .max_read_variant_count = _param.max_read_variant_count,
                               .max_read_variant_count_normalized = _param.max_read_variant_count_normalized,
                               .min_homopolymer_length = _param.min_homopolymer_length,
                               .sequencing_protocol = _param.sequencing_protocol,
                               .filter_homopolymer = _param.filter_homopolymer,
                               .tumor_sample_name = _param.tumor_sample_name,
                               .tumor_rg_ids = GetReadGroupIdsForSample(
                                   _alignment_reader_cache.Open(_param.bam_files, "r"), _param.tumor_sample_name),
                               .decode_yc = _param.decode_yc,
                               .min_base_type = _param.min_base_type,
                               .duplex_lowbq = _param.duplex_lowbq};
  _global_ctx.skip_variants_vcf = _param.skip_variants_vcf;

  // reference genome and region parameters
  _global_ctx.genome = _param.genome;
  _global_ctx.ref_seqs = &_ref_seqs;
  _global_ctx.bed_regions =
      _param.bed_file.has_value() ? GetChromIntervalMap(_param.bed_file).value() : ChromIntervalsMap{};
  _global_ctx.interest_regions =
      _param.interest_regions.has_value() ? _param.interest_regions.value() : ChromIntervalsMap{};

  _global_ctx.hdr = _hdr;
  std::tie(_global_ctx.vcf_info_metadata, _global_ctx.vcf_fmt_metadata) = _hdr->GetFieldMetadata();

  // Sex chromosome and normalization parameters — uses pre-computed median DP from the VCF pre-scan
  const auto& [sex, chr_x_name, chr_y_name, chr_x_par, chr_y_par, chr_med_dp] =
      GetSexInfo(_param, _pre_scan.chrom_median_dp);
  if (_param.normalize_features == FeatureNormalization::kMedianDp) {
    // Check whether the calculated median depth values are valid for normalization.
    ValidateChromMedianDepthForNormalization(chr_med_dp, FeatureNamesHaveSampleContext(_model_config.workflow));
  }
  // Store chromosome median DP when feature normalization is enabled or when the normalized DP ratio filter is active,
  // since the latter needs per-chromosome median DP values to compute the ratio.
  const bool needs_chrom_median_dp =
      _param.normalize_features == FeatureNormalization::kMedianDp || _param.min_dp_ratio > 0;
  _global_ctx.normalize_targets = needs_chrom_median_dp ? chr_med_dp : ChromMedianDepth{};

  // Pre-computed variant density from the VCF pre-scan (avoids boundary effects from per-region computation)
  _global_ctx.pre_scan_variants = &_pre_scan.chrom_variants;
  _global_ctx.pre_scan_density = &_pre_scan.chrom_variant_density;
  _global_ctx.normalize_scoring_features = (_param.normalize_features == FeatureNormalization::kMedianDp);
  _global_ctx.sex = sex;
  _global_ctx.chr_x_name = chr_x_name;
  _global_ctx.chr_y_name = chr_y_name;
  _global_ctx.chr_x_par = chr_x_par;
  _global_ctx.chr_y_par = chr_y_par;

  // Workflow specific parameters
  switch (_model_config.workflow) {
    case kGermlineMultiSample:
    case kGermline: {
      break;
    }
    case kTumorOnlyTe: {
      StrUnorderedSet block_list;
      if (_param.block_list) {
        block_list = LoadBlockList(*_param.block_list);
      }
      StrUnorderedSet hotspots;
      if (_param.hotspot_list) {
        hotspots = LoadHotspotVariants(*_param.hotspot_list);
      }
      _global_ctx.phased = _param.phased;
      _global_ctx.force_calls = GetChromIntervalMap(_param.forcecall_list);
      _global_ctx.hotspots = hotspots;
      _global_ctx.block_list = block_list;
      _global_ctx.min_allele_freq_threshold = _param.min_allele_freq_threshold;
      _global_ctx.weighted_counts_threshold = _param.weighted_counts_threshold;
      _global_ctx.hotspot_weighted_counts_threshold = _param.hotspot_weighted_counts_threshold;
      _global_ctx.ml_threshold = _param.ml_threshold;
      _global_ctx.hotspot_ml_threshold = _param.hotspot_ml_threshold;
      _global_ctx.min_phased_allele_freq = _param.min_phased_allele_freq;
      _global_ctx.max_phased_allele_freq = _param.max_phased_allele_freq;
      _global_ctx.min_alt_counts = _param.min_alt_counts;
      break;
    }
    case kGermlineTagging: {
      const auto tn_sample_idx = GetTumorNormalSampleIndexes(_hdr);
      if (tn_sample_idx.has_value()) {
        _global_ctx.vcf_tumor_index = tn_sample_idx->tumor_sample_idx;
        _global_ctx.vcf_normal_index = tn_sample_idx->normal_sample_idx;
      }
      _global_ctx.snv_min_ml_score = _param.snv_min_ml_score;
      _global_ctx.indel_min_ml_score = _param.indel_min_ml_score;
      _global_ctx.is_germline_tagging = true;
      break;
    }
    case kTumorNormalWgs: {
      _global_ctx.snv_min_ml_score = _param.snv_min_ml_score;
      _global_ctx.indel_min_ml_score = _param.indel_min_ml_score;
      _global_ctx.min_tumor_support = _param.min_tumor_support;
      _global_ctx.max_normal_support = _param.max_normal_support;
      _global_ctx.min_tumor_af = _param.min_tumor_af;
      _global_ctx.min_dp_ratio = _param.min_dp_ratio;
      _global_ctx.max_indel_size = _param.max_indel_size;
      const auto tn_sample_idx = GetTumorNormalSampleIndexes(_hdr);
      if (tn_sample_idx.has_value()) {
        _global_ctx.vcf_tumor_index = tn_sample_idx->tumor_sample_idx;
        _global_ctx.vcf_normal_index = tn_sample_idx->normal_sample_idx;
      }
      break;
    }
    default: {
      break;
    }
  }
}

/**
 * @brief Helper function to create score calculators based on workflow and model configuration.
 * @param param CLI parameters for `filter_variants`
 * @param model_config Model configuration for scoring
 * @return Vector of ScoreCalculator instances
 */
static vec<ScoreCalculator> CreateScoreCalculators(const FilterVariantsParam& param, const SVCConfig& model_config) {
  using enum Workflow;
  vec<ScoreCalculator> calculators;
  if (param.workflow == kGermline || param.workflow == kGermlineMultiSample) {
    calculators.emplace_back(param.snv_model,
                             model_config.snv_scoring_cols.size(),
                             model_config.snv_model_lgbm_prediction_params,
                             param.snv_shap_value_tsv.has_value());
    calculators.emplace_back(param.indel_model,
                             model_config.indel_scoring_cols.size(),
                             model_config.indel_model_lgbm_prediction_params,
                             param.indel_shap_value_tsv.has_value());
  } else {
    calculators.emplace_back(param.model,
                             model_config.scoring_cols.size(),
                             model_config.model_lgbm_prediction_params,
                             param.shap_value_tsv.has_value());
  }
  return calculators;
}

/**
 * @brief Run a single filter task for one genomic region.
 *
 * Dispatches to the appropriate workflow-specific filter method, collects the resulting
 * VCF records and SHAP rows, and pushes them to the writer queue. If the queue has been
 * shut down (e.g., due to a writer failure), the task returns immediately.
 *
 * @param worker Per-thread filter worker that owns the readers and calculators.
 * @param workflow Active workflow type (germline, tumor-only, tumor-normal, etc.).
 * @param region Genomic region to filter.
 * @param region_index Index of this region in the partitioned region list, used to
 *                     associate the result with the correct output file.
 * @param queue Thread-safe queue to push the completed result to.
 * @param progress Shared progress tracker for logging.
 * @throws error::Error if filtering fails for the region.
 */
static void RunFilterTask(FilterRegionClass& worker,
                          const Workflow workflow,
                          const TargetRegion& region,
                          const size_t region_index,
                          RegionOutputQueue& queue,
                          Progress& progress) {
  using enum Workflow;
  if (queue.IsShutdown()) {
    return;
  }
  try {
    progress.UpdateAndLog(log::LogLevel::kInfo);
    RegionResult result;
    switch (workflow) {
      case kGermlineMultiSample:
      case kGermline: {
        worker.FilterGermlineRegion(region, result);
        break;
      }
      case kGermlineTagging: {
        worker.FilterGermlineTaggingRegion(region, result);
        break;
      }
      case kTumorOnlyTe: {
        worker.FilterTumorOnlyTeRegion(region, result);
        break;
      }
      case kTumorNormalWgs: {
        worker.FilterTumorNormalRegion(region, result);
        break;
      }
      default: {
        break;
      }
    }
    queue.Push({region_index, std::move(result)});
  } catch (std::exception& e) {
    throw error::Error(
        "Error filtering variants in region '{}:{}-{}': {}", region.chrom, region.start, region.end, e.what());
  }
}

/**
 * @brief Writer thread body: consume completed regions and write each to a headerless BGZF file.
 *
 * Pops RegionOutput items from the queue, writes VCF records as raw BGZF data (no VCF header)
 * using vcf_format + bgzf_write, and appends any SHAP value rows to the corresponding TSV writers.
 * On failure, captures the exception and shuts down the queue to signal producer tasks to stop.
 *
 * @param num_regions Total number of regions expected. The function returns after consuming
 *                    this many items or when the queue is shut down.
 * @param queue Thread-safe queue to pop completed region results from.
 * @param region_files Ordered list of per-region output file paths, indexed by region_index.
 * @param hdr VCF header used to format records via VcfWriter::WriteBgzfRecord.
 * @param shap_writers Optional SHAP value TSV writers (combined, SNV-only, InDel-only).
 * @param exception_out Output parameter set to the caught exception if the writer fails.
 */
static void WriteRegionResults(const size_t num_regions,
                               RegionOutputQueue& queue,
                               const vec<fs::path>& region_files,
                               const io::VcfHeaderPtr& hdr,
                               const ShapWriters& shap_writers,
                               std::exception_ptr& exception_out) {
  try {
    size_t regions_written = 0;
    while (regions_written < num_regions) {
      auto item = queue.Pop();
      if (!item.has_value()) {
        break;
      }
      const auto& output = item.value();
      const auto& file_path = region_files.at(output.region_index);

      // Write this region's records as raw BGZF data (no VCF header).
      // vcf_format + bgzf_write bypasses bcf_write, which requires bcf_hdr_write first.
      if (!output.result.out_records.empty()) {
        const BgzfPtr bgzf(bgzf_open(file_path.c_str(), "w"));
        if (bgzf == nullptr) {
          throw error::Error("Failed to open BGZF file: {}", file_path.string());
        }
        for (const auto& record : output.result.out_records) {
          io::VcfWriter::WriteBgzfRecord(bgzf.get(), hdr, record);
        }
      }

      if (shap_writers.combined && !output.result.shap_value_rows.empty()) {
        shap_writers.combined->AppendRows(output.result.shap_value_rows);
      }
      if (shap_writers.snv && !output.result.snv_shap_value_rows.empty()) {
        shap_writers.snv->AppendRows(output.result.snv_shap_value_rows);
      }
      if (shap_writers.indel && !output.result.indel_shap_value_rows.empty()) {
        shap_writers.indel->AppendRows(output.result.indel_shap_value_rows);
      }
      ++regions_written;
    }
  } catch (const std::runtime_error&) {  // NOSONAR — std::runtime_error is the concrete type thrown by error::Error
    exception_out = std::current_exception();
    queue.Shutdown();
  }
}

void FilterVariantsClass::ParallelFiltering() {
  // Create a unique temp directory for per-region BGZF data files.
  // TmpDir uses mkdtemp for atomic, collision-free creation and installs signal handlers for cleanup.
  const TmpDir tmp_dir_guard(_param.vcf_output.parent_path(), ".svc_tmp");
  const auto& tmp_dir = tmp_dir_guard.Path();

  // Write the VCF header to a temp file. This triggers bcf_hdr_sync which rebuilds the
  // header's internal dictionary after AddFilterLine/AddFormatLine calls. Must happen
  // before SetGlobalContext() which reads the dictionary via GetFieldMetadata().
  const auto header_file = tmp_dir / "header.vcf.gz";
  {
    const io::VcfWriter header_writer(header_file, _hdr, io::NoIndex{});
    header_writer.WriteHeader();
  }

  // Initialize the global context — shared read-only state across all regions and tasks.
  SetGlobalContext();
  if (auto num_warnings = CheckVcfFeatureResources(
          _model_config.vcf_feature_cols, _global_ctx.interest_regions, _hdr, _param.pop_af_vcf.has_value());
      num_warnings > 0) {
    WarnAsErrorIfSet(
        "There were {} warnings while checking VCF feature resources; VCF feature extraction may not be accurate!",
        num_warnings);
  }

  const auto num_regions = _partitioned_regions.size();
  Progress progress(num_regions);

  // Build the ordered list of per-region temp file paths.
  vec<fs::path> region_files;
  region_files.reserve(num_regions);
  for (const auto idx : std::views::iota(size_t{0}, num_regions)) {
    region_files.emplace_back(tmp_dir / fmt::format("region_{}.vcf.gz", idx));
  }

  RegionOutputQueue queue;

  // Set up per-thread worker contexts.
  _workers.reserve(_param.threads);
  for (u32 i = 0; i < _param.threads; ++i) {
    try {
      auto worker_ctx = std::make_unique<WorkerContext>(
          _param.vcf_file, _param.pop_af_vcf, _param.bam_files, _alignment_reader_cache);
      worker_ctx->calculators = CreateScoreCalculators(_param, _model_config);
      _workers.emplace_back(_global_ctx, worker_ctx, _hdr);
    } catch (std::exception& e) {
      throw error::Error("Error creating WorkerContext: {}", e.what());
    }
  }

  // Build Taskflow graph: one filter task per region.
  tf::Taskflow taskflow;
  tf::Executor executor{_param.threads};
  const auto workflow = _model_config.workflow;

  for (size_t idx = 0; idx < num_regions; ++idx) {
    taskflow.emplace([&workers = _workers,
                      workflow,
                      &region = _partitioned_regions[idx],
                      region_idx = idx,
                      &queue,
                      &executor,
                      &progress]() {
      auto& worker = workers.at(static_cast<size_t>(executor.this_worker_id()));
      RunFilterTask(worker, workflow, region, region_idx, queue, progress);
    });
  }

  // Dedicated writer thread: pops completed regions and writes each to a headerless BGZF data file.
  const ShapWriters shap_writers{_shap_value_writer, _snv_shap_value_writer, _indel_shap_value_writer};
  std::exception_ptr writer_exception;
  std::jthread writer_thread([&num_regions, &queue, &region_files, &hdr = _hdr, &shap_writers, &writer_exception]() {
    WriteRegionResults(num_regions, queue, region_files, hdr, shap_writers, writer_exception);
  });

  Logging::Info("Filtering and writing '{}' regions...", progress.total_region_count);

  // Run all filter tasks. On failure, shut down the queue so the writer thread exits.
  std::exception_ptr filter_exception;
  try {
    executor.run(taskflow).get();
  } catch (const std::runtime_error&) {  // NOSONAR — tasks throw error::Error (std::runtime_error)
    filter_exception = std::current_exception();
    queue.Shutdown();
  }

  writer_thread.join();

  // Concatenate header + per-region BGZF data files at the block level (no record re-reading).
  if (!filter_exception && !writer_exception) {
    CreateParentDirectoryIfNotExists(_param.vcf_output);
    ConcatenateBgzfFiles(_param.vcf_output, header_file, region_files);
    if (bcf_index_build(_param.vcf_output.c_str(), 0) != 0) {
      Logging::Warn("Failed to build index for {}", _param.vcf_output.string());
    }
  }

  // Propagate exceptions — filter exceptions take priority.
  // tmp_dir is cleaned up by TmpDir destructor on scope exit.
  if (filter_exception) {
    std::rethrow_exception(filter_exception);
  }
  if (writer_exception) {
    std::rethrow_exception(writer_exception);
  }

  if (num_regions != progress.total_region_count) {
    throw error::Error("Number of regions {} and tasks executed {} are not the same!",
                       num_regions,
                       std::to_string(progress.total_region_count));
  }
}

void FilterVariantsClass::SetReferenceSequences() {
  StrSet chrom_names{};
  for (const auto& r : _partitioned_regions) {
    // Only store the chromosomes that have variants in the VCF file that need to be processed
    chrom_names.insert(r.chrom);
  }
  io::FastaReader fasta_reader(_param.genome);
  for (const auto& [chrom, length] : _hdr->GetContigLengths()) {
    if (chrom_names.contains(chrom)) {
      _ref_seqs[chrom] = fasta_reader.GetSequence(chrom, 0, static_cast<s32>(length));
    }
  }
}

void FilterVariantsClass::SetShapValueTsvWriters() {
  // Set up SHAP values output writers
  if (_param.shap_value_tsv.has_value()) {
    _shap_value_writer = std::make_shared<LockedTsvWriter>(_param.shap_value_tsv.value());
    WriteShapValueTsvHeader(*_shap_value_writer, _model_config.scoring_names, _param.command_line);
  }
  if (_param.snv_shap_value_tsv.has_value()) {
    _snv_shap_value_writer = std::make_shared<LockedTsvWriter>(_param.snv_shap_value_tsv.value());
    WriteShapValueTsvHeader(*_snv_shap_value_writer, _model_config.snv_scoring_names, _param.command_line);
  }
  if (_param.indel_shap_value_tsv.has_value()) {
    _indel_shap_value_writer = std::make_shared<LockedTsvWriter>(_param.indel_shap_value_tsv.value());
    WriteShapValueTsvHeader(*_indel_shap_value_writer, _model_config.indel_scoring_names, _param.command_line);
  }
}

void FilterVariantsClass::AddDuplexFormatLines() const {
  using enum io::FieldType;
  if (IsDuplexProtocol(_param.sequencing_protocol)) {
    _hdr->AddFormatLine({kDuplexConcordantCountsId, kDuplexConcordantCountsDesc, io::kNumberR, kInteger});
    _hdr->AddFormatLine({kDuplexSimplexCountsId, kDuplexSimplexCountsDesc, io::kNumberR, kInteger});
    _hdr->AddFormatLine({kDuplexDiscordantCountsId, kDuplexDiscordantCountsDesc, io::kNumberR, kInteger});
    _hdr->AddFormatLine({kDuplexLowbqCountsId, kDuplexLowbqCountsDesc, io::kNumberR, kInteger});
  }
}

void FilterVariantsClass::FilterGermline() {
  // Verify that the model files and scoring columns are compatible. The specified scoring columns must match the
  // features used in the model both in name and order. Any mismatch means the models are not compatible with the
  // requested workflow based on the provided config and will not produce reliable results.
  VerifyModelCompatibility(_param.snv_model, _model_config.snv_scoring_cols);
  VerifyModelCompatibility(_param.indel_model, _model_config.indel_scoring_cols);

  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  if (1 != reader.GetHeader()->GetNumSamples()) {
    // Germline filtering is designed for a single sample only.
    throw error::Error("VCF file must contain exactly one sample");
  }

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());
  // Add germline specific header information
  _hdr = reader.GetHeader()->Clone();
  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(io::GetCommandInfo(*_param.command_line));
  }
  _hdr->AddFilterLine({kFilteringPassId, kFilteringGermlinePassDesc});
  _hdr->AddFilterLine({kFilteringFailId, kFilteringGermlineFailDesc});
  _hdr->AddFilterLine({kFilteringMissingFeatureId, kFilteringMissingFeatureDesc});
  _hdr->AddFilterLine({kFilteringFalsePositiveId, kFilteringFalsePositiveDesc});
  _hdr->AddFilterLine({kFilteringMultialleleFormatId, kFilteringMultialleleFormatDesc});
  _hdr->AddFilterLine({kFilteringMultiallelePartnerId, kFilteringMultiallelePartnerDesc});
  _hdr->AddFilterLine({kFilteringMultialleleConflictId, kFilteringMultialleleConflictDesc});
  _hdr->AddFilterLine({kFilteringNonAcgtRefAltId, kFilteringNonAcgtRefAltDesc});
  _hdr->AddFormatLine({kGermlineMLId, kGermlineMLDesc, io::kNumberOne, io::FieldType::kInteger});
  _hdr->AddFormatLine({kMachineLearningId, kMachineLearningDesc, io::kNumberOne, io::FieldType::kFloat});
  AddDuplexFormatLines();
  _hdr->AddFormatLine({kGermlineGnomadAFId, kGermlineGnomadAFDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineRefAvgMapqId, kGermlineRefAvgMapqDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineAltAvgMapqId, kGermlineAltAvgMapqDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineRefAvgDistId, kGermlineRefAvgDistDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineAltAvgDistId, kGermlineAltAvgDistDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineDensity100BPId, kGermlineDensity100BPDesc, io::kNumberOne, io::FieldType::kInteger});

  // Single parallel VCF pre-scan: collects DP values, variant positions, and variant density
  _pre_scan = ParallelVcfPreScan(_param.vcf_file, _param.threads, std::nullopt);

  // Partition using DP-weighted regions for balanced parallel filtering
  _partitioned_regions = PartitionVcfRegionsByDp(_pre_scan, _param.threads);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  SetShapValueTsvWriters();

  // Filter the input VCF in parallel and write the output
  ParallelFiltering();

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::FilterGermlineTagging() {
  VerifyModelCompatibility(_param.model, _model_config.scoring_cols);

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());
  // Add germline specific header information
  _hdr = reader.GetHeader()->Clone();
  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(io::GetCommandInfo(*_param.command_line));
  }
  _hdr->AddFormatLine({kMachineLearningId, kMachineLearningDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine(
      {kGermlineTaggingInfoGermlineId, kGermlineTaggingInfoGermlineDesc, io::kNumberZero, io::FieldType::kFlag});
  _hdr->AddInfoLine(
      {kGermlineTaggingInfoSomaticId, kGermlineTaggingInfoSomaticDesc, io::kNumberZero, io::FieldType::kFlag});
  if (!_hdr->HasInfoField(kFilteringPassId)) {
    _hdr->AddFilterLine({kFilteringPassId, kFilteringPassDesc});
  }
  if (!_hdr->HasInfoField(kFilteringFailId)) {
    _hdr->AddFilterLine({kFilteringFailId, kFilteringFailDesc});
  }

  // Single parallel VCF pre-scan: collects DP values, variant positions, and variant density
  _pre_scan = ParallelVcfPreScan(_param.vcf_file, _param.threads, std::nullopt);

  // Partition using DP-weighted regions for balanced parallel filtering
  _partitioned_regions = PartitionVcfRegionsByDp(_pre_scan, _param.threads);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  SetShapValueTsvWriters();

  // Filter the input VCF in parallel and write the output
  ParallelFiltering();

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::FilterTumorNormal() {
  VerifyModelCompatibility(_param.model, _model_config.scoring_cols);

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());
  // Add germline specific header information
  _hdr = reader.GetHeader()->Clone();
  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(io::GetCommandInfo(*_param.command_line));
  }
  _hdr->AddFilterLine({kFilteringPassId, kFilteringGermlinePassDesc});
  _hdr->AddFilterLine({kFilteringFailId, kFilteringGermlineFailDesc});
  _hdr->AddFilterLine({kFilteringFailMinTumorAfId, kFilteringFailMinTumorAfDesc});
  _hdr->AddFilterLine({kFilteringMLScoreId, kFilteringMLScoreDesc});
  _hdr->AddFilterLine({kFilteringFailMinTumorSupportId, kFilteringFailMinTumorSupportDesc});
  _hdr->AddFilterLine({kFilteringFailMaxNormalSupportId, kFilteringFailMaxNormalSupportDesc});
  _hdr->AddFilterLine({kFilteringFailMinDpRatioId, kFilteringMinDpRatioDesc});
  _hdr->AddFilterLine({kFilteringFailMaxIndelSizeId, kFilteringFailMaxIndelSizeDesc});
  _hdr->AddFilterLine({kFilteringMissingFeatureId, kFilteringMissingFeatureDesc});
  _hdr->AddFilterLine({kFilteringNonAcgtRefAltId, kFilteringNonAcgtRefAltDesc});
  _hdr->AddInfoLine({kRefBQId, kRefBQDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine({kRefMQId, kRefMQDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine({kAltBQId, kAltBQDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine({kAltMQId, kAltMQDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine({kSubtypeId, kSubtypeDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine({kContextId, kContextDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine({kFieldPopaf, kPopafDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kMachineLearningId, kMachineLearningDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kBaseqQualId, kBaseQualDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kMapQualId, kMapQualDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kDistanceId, kDistanceDesc, io::kNumberOne, io::FieldType::kFloat});
  AddDuplexFormatLines();

  // Single parallel VCF pre-scan: collects DP values, variant positions, and variant density
  _pre_scan = ParallelVcfPreScan(_param.vcf_file, _param.threads, std::nullopt);

  // Partition using DP-weighted regions for balanced parallel filtering
  _partitioned_regions = PartitionVcfRegionsByDp(_pre_scan, _param.threads);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  SetShapValueTsvWriters();

  // Filter the input VCF in parallel and write the output
  ParallelFiltering();

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::FilterTumorOnlyTe() {
  VerifyModelCompatibility(_param.model, _model_config.scoring_cols);

  // load the VCF File
  Logging::Info("Loading VCF file {}", _param.vcf_file);
  io::VcfReader reader(_param.vcf_file);
  if (1 != reader.GetHeader()->GetNumSamples()) {
    throw error::Error("VCF file must contain exactly one sample");
  }

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());

  _hdr = reader.GetHeader()->Clone();
  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(io::GetCommandInfo(*_param.command_line));
  }
  // Append descriptions for the FILTER fields
  _hdr->AddFilterLine({kFilteringPassId, kFilteringPassDesc});
  _hdr->AddFilterLine({kFilteringMapQualityId, kFilteringMapQualityDesc});
  _hdr->AddFilterLine({kFilteringBlocklistedId, kFilteringBlocklistedDesc});
  _hdr->AddFilterLine({kFilteringMLScoreId, kFilteringMLScoreDesc});
  _hdr->AddFilterLine({kFilteringCountsId, kFilteringCountsDesc});
  _hdr->AddFilterLine({kFilteringBaseQualityId, kFilteringBaseQualityDesc});
  _hdr->AddFilterLine({kFilteringForcedId, kFilteringForcedDesc});
  _hdr->AddFilterLine({kFilteringAFId, kFilteringAFDesc});
  _hdr->AddFilterLine({kFilteringMinAltCountsId, kFilteringMinAltCountsDesc});

  // Append descriptions for the FORMAT fields
  _hdr->AddFormatLine({kWeightedCountsId, kWeightedCountsDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  _hdr->AddFormatLine({kMachineLearningId, kMachineLearningDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kNonDuplexCountsId, kNonDuplexCountsDesc, io::kNumberOne, io::FieldType::kInteger});
  _hdr->AddFormatLine({kDuplexCountsId, kDuplexCountsDesc, io::kNumberOne, io::FieldType::kInteger});
  _hdr->AddFormatLine({kStrandBiasMetricId, kStrandBiasMetricDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kPlusOnlyCountsId, kPlusOnlyCountsDesc, io::kNumberOne, io::FieldType::kInteger});
  _hdr->AddFormatLine({kMinusOnlyCountsId, kMinusOnlyCountsDesc, io::kNumberOne, io::FieldType::kInteger});
  _hdr->AddFormatLine({kSequenceContextId, kSequenceContextDesc, io::kNumberOne, io::FieldType::kString});
  _hdr->AddFormatLine({kPhysicalPhasedId, kPhysicalPhasedDesc, io::kNumberOne, io::FieldType::kString});
  _hdr->AddFormatLine({kHotspotId, kHotspotDesc, io::kNumberOne, io::FieldType::kInteger});
  AddDuplexFormatLines();

  // Load BED regions if provided
  const auto bed_regions = GetChromIntervalMap(_param.bed_file);

  // Pre-scan the full VCF (without BED filtering) so that median DP and variant density are
  // computed from all variants. BED filtering is applied only during region partitioning.
  _pre_scan = ParallelVcfPreScan(_param.vcf_file, _param.threads, std::nullopt);

  // Partition using DP-weighted regions, restricted to BED intervals when provided
  _partitioned_regions = PartitionVcfRegionsByDp(_pre_scan, _param.threads, 1, 1, bed_regions);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  SetShapValueTsvWriters();

  // Filter E2E by region and write the output
  ParallelFiltering();

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::CheckTumorNormalSampleName() const {
  if (!_param.tumor_sample_name.has_value() || _param.tumor_sample_name->empty()) {
    throw error::Error("Tumor sample name for read groups is not specified");
  }
  Logging::Info("Tumor sample name for read groups: {}", _param.tumor_sample_name.value());
}

void FilterVariantsClass::VerifyParameters() const {
  using enum Workflow;
  switch (_model_config.workflow) {
    case kGermlineTagging:
    case kTumorNormalWgs: {
      CheckTumorNormalSampleName();
      break;
    }
    case kTumorOnlyTe:
    case kGermline:
    case kGermlineMultiSample:
      break;
    default:
      throw error::Error("Unsupported workflow");
  }

  if (_param.bam_files.empty()) {
    throw error::Error("No BAM file provided");
  }
}

/**
 * The main entry point for the filter_variants module. This function performs checks to ensure the correct number of
 * models are passed based on the specified workflow before calling a workflow specific filtering function.
 * @param param A FilterVariantsParam struct that contains required input/output filepaths and filtering settings
 */
void FilterVariants(const FilterVariantsParam& param) {
  using enum Workflow;
  // This function acts as the entry point for the filter-variants tool. Based on the workflow the respective filtering
  // process is run

  // Create Filtering object
  FilterVariantsClass filtering(param);

  // Verify that the correct number of models have been passed based on the workflow, and that BAM input has been passed
  filtering.VerifyParameters();

  switch (param.workflow) {
    case kGermlineMultiSample:
    case kGermline: {
      filtering.FilterGermline();
      break;
    }
    case kGermlineTagging: {
      filtering.FilterGermlineTagging();
      break;
    }
    case kTumorOnlyTe: {
      filtering.FilterTumorOnlyTe();
      break;
    }
    case kTumorNormalWgs: {
      filtering.FilterTumorNormal();
      break;
    }
    default:
      throw error::Error("Specified Workflow has not yet been implemented");
  }
  Logging::Info("Filtered all variants");
}

}  // namespace xoos::svc
