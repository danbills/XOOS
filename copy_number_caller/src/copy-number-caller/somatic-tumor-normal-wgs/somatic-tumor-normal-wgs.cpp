#include "copy-number-caller/somatic-tumor-normal-wgs/somatic-tumor-normal-wgs.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>

#include <taskflow/algorithm/for_each.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/util/file-functions.h>

#include "baits.h"
#include "copy-number-caller/common/bam-to-corrected-coverage.h"
#include "copy-number-caller/common/check-chr-y-par.h"
#include "copy-number-caller/common/coverage-check.h"
#include "copy-number-caller/common/write-empty-outputs.h"
#include "coverage.h"
#include "filter-observations.h"
#include "io/baits-io.h"
#include "io/copy-number-caller-default-filenames.h"
#include "io/write-baf-values.h"
#include "io/write-igv-xml.h"
#include "io/write-logr.h"
#include "io/write-metrics.h"
#include "io/write-segments.h"
#include "likelihood/likelihood-model.h"
#include "likelihood/likelihood.h"
#include "mapq-utils.h"
#include "misc/sample-metadata-options.h"
#include "observations.h"
#include "purity-ploidy-search/purity-ploidy-search.h"
#include "seg-to-vcf/seg-to-vcf.h"
#include "segmentation/parent-specific-binary-segmentation-taskflow-graph.h"
#include "segmentation/read-segments.h"
#include "two-sample-logr/tumor-normal.h"
#include "vcf-purity-source.h"

namespace xoos::cnc {

constexpr u8 kDefaultPloidyForPredictTumorPurityFromSomaticVariants = 2;
constexpr f64 kDefaultMinTumorPurity = 0;
constexpr f64 kDefaultMaxTumorPurity = 0.99;
constexpr f64 kTumorPurityClampEpsilon = 1e-12;

struct SomaticTumorNormalWGSData {
  Observations logrs;
  Observations mapqs;
  Sex sex = Sex::kUnknown;
  Observations ref_depths;
  Observations alt_depths;
  std::vector<GenomicSegment> segments;
  double tumor_purity = 0.99;
  double tumor_ploidy = 2;
  VcfPuritySource vcf_purity_source = VcfPuritySource::kNone;
  std::optional<TotalCopyNumberPrior> total_copy_number_prior;
  std::optional<CoverageCheckResult> coverage_check_result;
};

SomaticTumorNormalWGSOutputPaths SetupDefaultSomaticTumorNormalWGSOutputPaths(const CopyNumberCallerOptions& options) {
  SomaticTumorNormalWGSOutputPaths paths{
      .likelihood_out = options.output_dir / kDefaultSomaticCNCallsetSegOutput,
      .vcf_out = options.output_dir / kDefaultSomaticCNCallsetVcfOutput,
      .metrics_out = options.output_dir / kDefaultMetricsOutput,
      .igv_xml_out = options.output_dir / kDefaultIgvXmlOutput,
      .logrs_out = options.output_dir / kDefaultLogRsOutput,
      .logrs_bw_out = options.output_dir / kDefaultLogRsBwOutput,
      .logr_segments_out = options.output_dir / kDefaultLogRsSegOutput,
      .purity_ploidy_grid_out = options.output_dir / kDefaultPurityPloidyGridOutput,
      .mapping_qualities_out = options.output_dir / kDefaultMapQOutput,
      .augmented_baits_out = options.output_dir / kDefaultAugmentedBaitsOutput,
      .tumor_coverage_out = options.output_dir / kDefaultTumorCoverageOutput,
      .tumor_corrected_coverage_out = options.output_dir / kDefaultTumorCorrectedCoverageOutput,
      .normal_coverage_out = options.output_dir / kDefaultNormalCoverageOutput,
      .normal_corrected_coverage_out = options.output_dir / kDefaultNormalCorrectedCoverageOutput,
      .baf_out = options.output_dir / kDefaultBafOutput,
      .baf_bw_out = options.output_dir / kDefaultBafBwOutput,
  };
  if (options.save_baf_segments) {
    paths.baf_segments_out = options.output_dir / kDefaultBafSegOutput;
  }
  if (options.save_taskflow_graph) {
    paths.taskflow_graph_out = options.output_dir / kDefaultTaskflowGraphOutput;
  }
  std::vector<fs::path> output_paths = {
      paths.likelihood_out,
      paths.vcf_out,
      paths.metrics_out,
      paths.igv_xml_out,
      paths.logrs_out,
      paths.logrs_bw_out,
      paths.logr_segments_out,
      paths.purity_ploidy_grid_out,
      paths.mapping_qualities_out,
      paths.augmented_baits_out,
      paths.tumor_coverage_out,
      paths.tumor_corrected_coverage_out,
      paths.normal_coverage_out,
      paths.normal_corrected_coverage_out,
      paths.baf_out,
      paths.baf_bw_out,
  };
  if (paths.baf_segments_out.has_value()) {
    output_paths.push_back(paths.baf_segments_out.value());
  }
  if (paths.taskflow_graph_out.has_value()) {
    output_paths.push_back(paths.taskflow_graph_out.value());
  }
  file::CheckFilePermissionsAndOutputPathExistence({}, output_paths);
  return paths;
}

/*
 * @brief Taskflow subroutine to calculate coverage and logR values
 * @param subflow Subflow
 * @param options Copy number caller options
 * @param baits Bait records for coverage calculation
 * @param data Data structure to hold results
 */
static void CoverageAndLogRsTask(tf::Subflow& subflow,
                                 const CopyNumberCallerOptions& options,
                                 const SomaticTumorNormalWGSOutputPaths& paths,
                                 const BaitRecords& baits,
                                 SomaticTumorNormalWGSData& data) {
  Logging::Info("calculating and correcting coverage for {} and {}",
                options.tumor_bam_fname.string(),
                options.normal_bam_fname.string());
  // tumor coverage calculation
  BamToCorrectedCoverageTaskFlowGraph tumor_bam_corrected_coverage_tf_graph(
      options.tumor_bam_fname,
      baits,
      paths.tumor_coverage_out,
      paths.tumor_corrected_coverage_out,
      options.calculate_coverage_options.exclude_flags,
      options.calculate_coverage_options.ignore_DN,
      options.gc_correct_options.first_span,
      options.command_line_info);
  tf::Task tumor_coverage_task =
      subflow.composed_of(tumor_bam_corrected_coverage_tf_graph).name("tumor coverage and correction");
  // normal coverage calculation
  BamToCorrectedCoverageTaskFlowGraph normal_corrected_coverage_tf_graph(
      options.normal_bam_fname,
      baits,
      paths.normal_coverage_out,
      paths.normal_corrected_coverage_out,
      options.calculate_coverage_options.exclude_flags,
      options.calculate_coverage_options.ignore_DN,
      options.gc_correct_options.first_span,
      options.command_line_info);
  tf::Task normal_coverage_task =
      subflow.composed_of(normal_corrected_coverage_tf_graph).name("normal coverage and correction");
  // we have to make sure that both coverage/correction tasks are performed in succession, in order to avoid thread
  // thrashing.
  tumor_coverage_task.succeed(normal_coverage_task);
  subflow.join();
  // Check for low coverage early exit on either tumor or normal
  if (tumor_bam_corrected_coverage_tf_graph.IsLowCoverage()) {
    data.coverage_check_result = tumor_bam_corrected_coverage_tf_graph.GetCoverageCheckResult();
    return;
  }
  if (normal_corrected_coverage_tf_graph.IsLowCoverage()) {
    data.coverage_check_result = normal_corrected_coverage_tf_graph.GetCoverageCheckResult();
    return;
  }
  // gather data and calculate log ratios
  auto normal_coverage = normal_corrected_coverage_tf_graph.GetResult();
  auto tumor_coverage = tumor_bam_corrected_coverage_tf_graph.GetResult();
  data.logrs =
      ProcessTumorNormal(tumor_coverage, normal_coverage, options.normal_min_coverage, options.tumor_min_coverage);
  data.sex = options.sample_metadata_options.sex.has_value() ? options.sample_metadata_options.sex.value()
                                                             : normal_coverage.PredictSex();
  // populate return values and write to disk if necessary
  data.mapqs = GetAvgMapqsFromCoverage(tumor_coverage);
  WriteLogRFiles(
      data.logrs, paths.logrs_out, paths.logrs_bw_out, options.reference_genome_fai_fname, options.command_line_info);
  std::ofstream mapq_ofs(paths.mapping_qualities_out);
  data.mapqs.Write(mapq_ofs, kMeanMapqColumnName, true, options.command_line_info);
}

/*
 * @brief Estimate tumor purity from somatic SNV VAFs and update vcf_purity_source on success
 * @param vcf_purity_source purity source indicator
 * @param somatic_vafs non-empty vector of somatic SNV VAFs to use for tumor purity estimation
 * @param tumor_purity estimated purity
 */
static void PredictTumorPurityFromSomaticVariants(VcfPuritySource& vcf_purity_source,
                                                  const vec<f64>& somatic_vafs,
                                                  f64& tumor_purity) {
  // calculate median VAF
  std::vector<f64> sorted_vafs(somatic_vafs.begin(), somatic_vafs.end());
  const auto median_idx = sorted_vafs.size() / 2;
  const auto median_it = std::next(sorted_vafs.begin(), static_cast<std::ptrdiff_t>(median_idx));
  std::ranges::nth_element(sorted_vafs, median_it);
  const f64 median_vaf = *median_it;
  // predict tumor purity based on median VAF (and clamp to [0, 1])
  const f64 unclamped_tumor_purity = kDefaultPloidyForPredictTumorPurityFromSomaticVariants * median_vaf;
  const f64 clamped_tumor_purity = std::clamp(unclamped_tumor_purity, kDefaultMinTumorPurity, kDefaultMaxTumorPurity);
  if (std::fabs(clamped_tumor_purity - unclamped_tumor_purity) > kTumorPurityClampEpsilon) {
    Logging::Warn(
        "Predicted tumor purity from somatic SNV VAFs was out of bounds and has been clamped: "
        "raw purity = {:.2f}, clamped purity = {:.2f}.",
        unclamped_tumor_purity,
        clamped_tumor_purity);
  }
  tumor_purity = clamped_tumor_purity;
  Logging::Info("Reestimated tumor purity from somatic SNV VAFs: median VAF = {:.2f}, predicted purity = {:.2f}.",
                median_vaf,
                tumor_purity);
  vcf_purity_source = VcfPuritySource::kSomaticSnv;
}

/*
 * @brief Taskflow subroutine to segment and call copy number
 * @param subflow Subflow
 * @param options Copy number caller options
 * @param baits Bait records for coverage calculation
 * @param seed_segments Seed segments for segmentation
 * @param ref_depths Reference allele depths from VCF
 * @param alt_depths Alternate allele depths from VCF
 * @param data Data structure to hold results
 */
static void SegmentationAndCallingTask(tf::Subflow& subflow,
                                       const CopyNumberCallerOptions& options,
                                       const SomaticTumorNormalWGSOutputPaths& paths,
                                       const BaitRecords& baits,
                                       const std::vector<GenomicSegment>& seed_segments,
                                       const RefAltObservations& ref_alt_depths,
                                       SomaticTumorNormalWGSData& data) {
  SegmentationOptions segmentation_options = options.segmentation_options;
  SampleMetadataOptions sample_metadata_options = options.sample_metadata_options;
  sample_metadata_options.sex = data.sex;
  // filter the SNPs to only those that overlap LogR intervals. Non-intersecting SNPs are irrelevant to the algorithm
  data.ref_depths = GetOverlappingObservations(data.logrs, ref_alt_depths.ref_obvs);
  data.alt_depths = GetOverlappingObservations(data.logrs, ref_alt_depths.alt_obvs);
  WriteBafFiles(data.ref_depths,
                data.alt_depths,
                paths.baf_out,
                paths.baf_bw_out,
                options.reference_genome_fai_fname,
                options.command_line_info);
  // calculate DH values (decrease in heterozygosity) from depths. These will be used for 2nd round of segmentation
  auto filtered_dh_vals = GetDhFromDepths(data.ref_depths, data.alt_depths);
  segmentation::ParentSpecificBinarySegmentationTaskflowGraph segment_graph(
      data.logrs, filtered_dh_vals, seed_segments, segmentation_options, sample_metadata_options);
  subflow.composed_of(segment_graph);
  subflow.join();
  WriteLogRSegments(segment_graph.GetLogRSegments(),
                    paths.logr_segments_out,
                    data.logrs,
                    options.sample_metadata_options.sample_id,
                    data.sex,
                    options.command_line_info);
  data.segments = segment_graph.GetResult();
  for (auto& seg : data.segments) {
    seg.id = options.sample_metadata_options.sample_id;
    seg.sex = data.sex;
  }
  if (paths.baf_segments_out.has_value()) {
    WriteBafSegments(data.segments, paths.baf_segments_out.value(), options.command_line_info);
  }
  Logging::Info("found {} segments", data.segments.size());
  // calculate tumor purity/ploidy if user asks for the grid
  if (sample_metadata_options.purity.has_value() && sample_metadata_options.ploidy.has_value()) {
    data.tumor_purity = sample_metadata_options.purity.value();
    data.tumor_ploidy = sample_metadata_options.ploidy.value();
    data.vcf_purity_source = VcfPuritySource::kUserInput;
  } else {
    const auto purity_ploidy_out = PurityPloidySearchWithWGDQC(data.segments,
                                                               data.logrs,
                                                               data.ref_depths,
                                                               data.alt_depths,
                                                               options.purity_ploidy_search_options,
                                                               paths.purity_ploidy_grid_out,
                                                               data.sex,
                                                               options.command_line_info);
    data.tumor_purity = purity_ploidy_out.purity;
    data.tumor_ploidy = purity_ploidy_out.ploidy;
    data.vcf_purity_source = VcfPuritySource::kScna;
    data.total_copy_number_prior = purity_ploidy_out.total_copy_number_prior;
    if (data.tumor_purity < kLowPredictedPurity) {
      Logging::Info(
          "Estimated or provided tumor purity is quite low ({:.2f}). Copy number calls may be unreliable at low purity "
          "levels.",
          data.tumor_purity);
      if (!ref_alt_depths.somatic_vafs.has_value() || ref_alt_depths.somatic_vafs->empty()) {
        Logging::Info(
            "Consider using a VCF with somatic SNV calls for better tumor purity estimation. No somatic SNVs found in "
            "VCF, skipping tumor purity prediction from somatic variants.");
      } else {
        PredictTumorPurityFromSomaticVariants(
            data.vcf_purity_source, ref_alt_depths.somatic_vafs.value(), data.tumor_purity);
        data.tumor_ploidy = kDefaultPloidyForPredictTumorPurityFromSomaticVariants;
      }
    }
  }
  // call copy numbers
  LikelihoodModel likelihood_model = LikelihoodModel::kSerialSummarized;
  LikelihoodOptions likelihood_options = options.likelihood_options;
  sample_metadata_options.purity = data.tumor_purity;
  sample_metadata_options.ploidy = data.tumor_ploidy;
  Logging::Info("Calculating likelihood and calling copy numbers using tumor purity {:.2f} and ploidy {:.2f}",
                data.tumor_purity,
                data.tumor_ploidy);
  data.segments = CalculateLikelihoods(data.segments,
                                       data.logrs,
                                       data.ref_depths,
                                       data.alt_depths,
                                       data.mapqs,
                                       likelihood_model,
                                       likelihood_options,
                                       sample_metadata_options);
}

CopyNumberCallerRet SomaticTumorNormalWGS(const BaitRecords& baits,
                                          const std::vector<GenomicSegment>& seed_segments,
                                          const RefAltObservations& ref_alt_depths,
                                          const CopyNumberCallerOptions& options,
                                          const SomaticTumorNormalWGSOutputPaths& paths) {
  // check conflicting options
  const bool has_purity = options.sample_metadata_options.purity.has_value();
  const bool has_ploidy = options.sample_metadata_options.ploidy.has_value();
  // 1. Check for Partial Parameters (e.g., Purity but no Ploidy)
  if (has_purity != has_ploidy) {
    Logging::Error("Must specify BOTH --tumor-purity and --tumor-ploidy, or NEITHER.");
    throw std::runtime_error("invalid options");
  }
  SomaticTumorNormalWGSData data;
  tf::Taskflow taskflow;
  tf::Task coverage_and_logrs_task = taskflow
                                         .emplace([&baits, &options, &paths, &data](tf::Subflow& subflow) {
                                           CoverageAndLogRsTask(subflow, options, paths, baits, data);
                                         })
                                         .name("coverage and logrs");
  tf::Task segmentation_and_likelihood_task =
      taskflow.emplace([&options, &paths, &baits, &seed_segments, &ref_alt_depths, &data](tf::Subflow& subflow) {
        if (IsLowCoverage(data.coverage_check_result)) {
          return;
        }
        SegmentationAndCallingTask(subflow, options, paths, baits, seed_segments, ref_alt_depths, data);
      });
  coverage_and_logrs_task.precede(segmentation_and_likelihood_task);
  if (paths.taskflow_graph_out.has_value()) {
    std::ofstream ofs(paths.taskflow_graph_out.value());
    taskflow.dump(ofs);
  }
  Logging::Info("using {} threads", options.threads);
  tf::Executor executor(options.threads);
  executor.run(taskflow).get();

  // Handle low coverage early exit
  if (IsLowCoverage(data.coverage_check_result)) {
    CopyNumberCallerRet ret{};
    ret.coverage_check_result = data.coverage_check_result;
    return ret;
  }

  // add some metadata to all the segments
  for (auto& seg : data.segments) {
    seg.id = options.sample_metadata_options.sample_id;
    seg.sex = data.sex;
  }
  return {data.ref_depths,
          data.alt_depths,
          data.segments,
          data.sex,
          data.tumor_purity,
          data.tumor_ploidy,
          data.vcf_purity_source,
          data.total_copy_number_prior};
}

void SomaticTumorNormalWGSMain(const CopyNumberCallerOptions& options) {
  if (options.mode != CopyNumberCallerModes::kSomaticTumorNormalWGS) {
    Logging::Error("invalid CopyNumberCaller mode for SomaticTumorNormalWGSMain");
    throw std::runtime_error("invalid CopyNumberCaller mode for SomaticTumorNormalWGSMain");
  }
  const auto paths = SetupDefaultSomaticTumorNormalWGSOutputPaths(options);
  BaitRecords baits = LoadOrGenerateWholeGenomeIntervals(options, paths.augmented_baits_out);
  io::VcfReader vcf_reader_for_baf(options.vcf_fname.value());
  const auto ref_alt_depths =
      GetDepthsFromVcf(vcf_reader_for_baf, false, true, options.baf_filter_options, options.sample_metadata_options);

  // Graceful exit: no variants remain after depth/BAF filtering
  if (ref_alt_depths.ref_obvs.contigs.empty() || ref_alt_depths.alt_obvs.contigs.empty()) {
    Logging::Warn("No variants found in VCF after applying filters. Writing empty outputs.");
    WriteSomaticEmptyOutputs(paths, std::nullopt, baits, options, ref_alt_depths.vcf_check);
    return;
  }

  // Graceful exit: low het-SNP fraction in somatic-flagged VCF
  if (IsVcfInsufficient(ref_alt_depths.vcf_check)) {
    const auto& vc = *ref_alt_depths.vcf_check;
    Logging::Warn(
        "The input VCF was detected as a somatic-flagged VCF (FILTER == PASS indicates a somatic variant), but "
        "only {:.1f}% of its evaluated variants are heterozygous germline variants, which is below the expected "
        "threshold of {:.0f}%. Writing empty outputs. To force analysis, use "
        "`--force-enable-somatic-variant-parsing`.",
        vc.het_snp_fraction.value_or(0.0) * 100.0,
        kDefaultMinHetSnpFraction * 100.0);
    WriteSomaticEmptyOutputs(paths, std::nullopt, baits, options, ref_alt_depths.vcf_check);
    return;
  }
  const std::vector<GenomicSegment> seed_segments =
      segmentation::ReadSegments(options.seed_segments_fname, SegmentType::kSeed);
  baits = RemovePARIfChrYPARUnmasked(options.normal_bam_fname, baits, seed_segments);
  baits = RemovePARIfChrYPARUnmasked(options.tumor_bam_fname, baits, seed_segments);
  const auto ret = SomaticTumorNormalWGS(baits, seed_segments, ref_alt_depths, options, paths);

  // Handle low coverage early exit — write header-only outputs and return
  if (IsLowCoverage(ret.coverage_check_result)) {
    WriteSomaticEmptyOutputs(paths, *ret.coverage_check_result, baits, options);
    return;
  }

  WriteSegments(
      paths.likelihood_out,
      ret.segments,
      options.vcf_fname.has_value() ? SegmentType::kSomaticWithBafLikelihood : SegmentType::kSomaticNoBafLikelihood,
      options.command_line_info);
  bool has_high_extreme_baf_proportion = false;
  if (ret.tumor_purity.has_value()) {
    has_high_extreme_baf_proportion =
        CheckExtremeBafProportionAndWarn(ret.ref_depths, ret.alt_depths, ret.tumor_purity.value());
  }
  segmentation::WriteSegmentsToVcf(paths.vcf_out,
                                   ret.segments,
                                   {.seq_lengths = baits.GetSeqLengths(),
                                    .seq_order = GetContigOrder(options.reference_genome_fai_fname),
                                    .sample_id = options.sample_metadata_options.sample_id,
                                    .reference_file = baits.GetReferenceFile(),
                                    .sex = ret.sex,
                                    .mode = options.mode,
                                    .purity = ret.tumor_purity,
                                    .ploidy = ret.tumor_ploidy,
                                    .vcf_purity_source = ret.vcf_purity_source,
                                    .total_copy_number_prior = ret.total_copy_number_prior,
                                    .command_line_info = options.command_line_info,
                                    .has_high_extreme_baf_proportion = has_high_extreme_baf_proportion});
  WriteMetricsFile(paths.metrics_out,
                   baits,
                   ret.ref_depths.contigs.size(),
                   ret.sex,
                   ret.segments,
                   ret.tumor_purity,
                   ret.tumor_ploidy,
                   options.mode,
                   options.command_line_info,
                   std::nullopt);
  WriteIGVXML(paths.igv_xml_out, std::make_optional(paths.logrs_bw_out), std::make_optional(paths.baf_bw_out));
}

}  // namespace xoos::cnc
