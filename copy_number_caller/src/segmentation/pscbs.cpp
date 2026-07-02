#include "segmentation/pscbs.h"

#include <fstream>

#include <taskflow/algorithm/for_each.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>
#include <xoos/util/file-functions.h>

#include "copy-number-caller/copy-number-caller-options.h"
#include "filter-observations.h"
#include "io/copy-number-caller-default-filenames.h"
#include "io/write-baf-values.h"
#include "io/write-bigwig.h"
#include "io/write-segments.h"
#include "misc/sample-metadata-options.h"
#include "observations.h"
#include "segmentation/parent-specific-binary-segmentation-taskflow-graph.h"
#include "segmentation/read-segments.h"
#include "segmentation/segment-type.h"
#include "segmentation/segmentation-options.h"
#include "sex.h"

namespace xoos::cnc::segmentation {

static SegmentType GetSegmentType(const CopyNumberCallerOptions& options) {
  if (options.vcf_fname.has_value()) {
    return SegmentType::kBaf;
  } else {
    return SegmentType::kLogROnly;
  }
}

/**
 * @brief entry-function for the Parent-Specific-Binary-Segmentation algorithm
 * @param cli_opts contains file names etc to be parsed for input to ParentSpecificBinarySegmentation function
 */
void ParentSpecificBinarySegmentationMain(const CopyNumberCallerOptions& options) {
  const auto logr_segments_out = options.output_dir / kDefaultLogRsSegOutput;
  const auto seg_out = options.output_dir / kDefaultSegmentationSegOutput;
  std::vector<fs::path> output_paths = {logr_segments_out, seg_out};
  if (options.vcf_fname.has_value()) {
    output_paths.push_back(options.output_dir / kDefaultBafOutput);
    if (!options.reference_genome_fai_fname.empty()) {
      output_paths.push_back(options.output_dir / kDefaultBafBwOutput);
    }
  }
  file::CheckFilePermissionsAndOutputPathExistence({}, output_paths);
  std::ifstream logrs_stream(options.logrs_fname.value());
  Observations logrs = ReadObservations(logrs_stream);
  Observations dh_values;
  if (options.vcf_fname.has_value()) {
    io::VcfReader vcf_reader(options.vcf_fname.value().string());
    Observations ref_depths;
    Observations alt_depths;
    if (options.sample_metadata_options.normal_sample_name.has_value() &&
        options.sample_metadata_options.tumor_sample_name.has_value()) {
      auto ref_alts =
          GetDepthsFromVcf(vcf_reader, false, true, options.baf_filter_options, options.sample_metadata_options);
      ref_depths = std::move(ref_alts.ref_obvs);
      alt_depths = std::move(ref_alts.alt_obvs);
    } else if (options.sample_metadata_options.normal_sample_name.has_value()) {
      auto ref_alts =
          GetDepthsFromVcf(vcf_reader, false, false, options.baf_filter_options, options.sample_metadata_options);
      ref_depths = std::move(ref_alts.ref_obvs);
      alt_depths = std::move(ref_alts.alt_obvs);
    } else if (options.sample_metadata_options.tumor_sample_name.has_value()) {
      SampleMetadataOptions tumor_as_normal = options.sample_metadata_options;
      tumor_as_normal.normal_sample_name = options.sample_metadata_options.tumor_sample_name;
      BafFilterOptions tumor_filter{options.baf_filter_options};
      tumor_filter.normal_sample_min_depth = options.baf_filter_options.tumor_sample_min_depth;
      tumor_filter.normal_sample_min_baf = options.baf_filter_options.tumor_sample_min_baf;
      tumor_filter.normal_sample_max_baf = options.baf_filter_options.tumor_sample_max_baf;
      auto ref_alts = GetDepthsFromVcf(vcf_reader, false, false, tumor_filter, tumor_as_normal);
      ref_depths = std::move(ref_alts.ref_obvs);
      alt_depths = std::move(ref_alts.alt_obvs);
    }
    ref_depths = GetOverlappingObservations(logrs, ref_depths);
    alt_depths = GetOverlappingObservations(logrs, alt_depths);
    dh_values = GetDhFromDepths(ref_depths, alt_depths);
    dh_values.SetBAFSegObvsStatus(true);
    const auto baf_out = options.output_dir / kDefaultBafOutput;
    Logging::Info("Writing BAF values to {}", baf_out.string());
    WriteBafValues(ref_depths, alt_depths, baf_out, options.command_line_info);
    if (!options.reference_genome_fai_fname.empty()) {
      const auto baf_bw_out = options.output_dir / kDefaultBafBwOutput;
      Logging::Info("Writing BAF values to {}", baf_bw_out.string());
      WriteBafValuesToBigWig(ref_depths, alt_depths, baf_bw_out, options.reference_genome_fai_fname);
    } else {
      Logging::Info("Reference FAI file was not provided and is required to write BAF values to BigWig");
    }
  } else if (options.ad_fname.has_value()) {
    const auto ref_alts = GetDepthsFromFile(options.ad_fname.value());
    dh_values = GetDhFromDepths(ref_alts.ref_obvs, ref_alts.alt_obvs);
    dh_values = GetOverlappingObservations(logrs, dh_values);
    dh_values.SetBAFSegObvsStatus(true);
  }
  std::vector<GenomicSegment> seed_segments;
  std::ifstream seed_segments_stream;
  const auto& command_line_info = options.command_line_info;
  if (!options.seed_segments_fname.empty()) {
    seed_segments = ReadSegments(options.seed_segments_fname, SegmentType::kSeed);
  } else {
    Logging::Error("Seed segments are required input!");
    throw std::runtime_error("missing seed segments");
  }
  std::vector<GenomicSegment> segments = ParentSpecificBinarySegmentation(logrs,
                                                                          dh_values,
                                                                          seed_segments,
                                                                          options.segmentation_options,
                                                                          options.sample_metadata_options,
                                                                          options.threads,
                                                                          logr_segments_out,
                                                                          command_line_info);
  for (auto& seg : segments) {
    seg.id = options.sample_metadata_options.sample_id;
    seg.sex = options.sample_metadata_options.sex;
  }
  SegmentType mode = GetSegmentType(options);
  WriteSegments(seg_out, segments, mode, command_line_info);
}

/**
 * @brief implementation of the Parent-Specific-Binary-Segmentation algorithm
 * described in Olshen et. al 2011
 * (https://doi.org/10.1093/bioinformatics/btr329). Briefly, this algorithm
 * performs Circular Binary Segmentation (CBS) first on the tumor-normal logr
 * values, then feeds the resulting segments into a second CBS step based on
 * B-allele fractions at germline heterozygous variants
 * @param logrs tumor-normal(panel) log-ratios for the sample
 * @param dh_vals (optional) the dh (inverted BAF) values for the target
 * @param seed_segments  seed segments (if empty, will not be used). Common seeds are partitions dividing chromosomes
 by centromeres/telomeres
 * @param max_p maximum p-value for calling a segment significant (default: 0.01)
 * @param n_permutations  number of permutations for significance testing
 * @param min_obs_per_segment minimum number of observations for a segment
 * @param n_threads number of threads to use
 */
std::vector<GenomicSegment> ParentSpecificBinarySegmentation(const Observations& logrs,
                                                             const Observations& dh_vals,
                                                             const std::vector<GenomicSegment>& seed_segments,
                                                             const SegmentationOptions& segmentation_options,
                                                             const SampleMetadataOptions& sample_metadata_options,
                                                             const size_t n_threads,
                                                             const fs::path& logr_segments_out,
                                                             const io::CommandLineInfo& command_line_info) {
  tf::Taskflow taskflow;
  ParentSpecificBinarySegmentationTaskflowGraph pscbs_graph(
      logrs, dh_vals, seed_segments, segmentation_options, sample_metadata_options);
  tf::Task pscbs_task = taskflow.composed_of(pscbs_graph);
  std::vector<GenomicSegment> res;
  tf::Task pscbs_res_task = taskflow.emplace(
      [&res, &logrs, &dh_vals, &pscbs_graph, &logr_segments_out, &command_line_info, &sample_metadata_options]() {
        res = pscbs_graph.GetResult();
        // repopulate segments with original array data
        PopulateGenomicSegmentOptionalFields(res, logrs);
        if (!dh_vals.starts.empty()) {
          PopulateGenomicSegmentOptionalFields(res, dh_vals);
        }
        if (!logr_segments_out.empty()) {
          WriteLogRSegments(pscbs_graph.GetLogRSegments(),
                            logr_segments_out,
                            logrs,
                            sample_metadata_options.sample_id,
                            sample_metadata_options.sex.value_or(Sex::kUnknown),
                            command_line_info);
        }
      });
  pscbs_task.precede(pscbs_res_task);
  tf::Executor executor(n_threads);
  executor.run(taskflow).get();
  return res;
}

}  // namespace xoos::cnc::segmentation
