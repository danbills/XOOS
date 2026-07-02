#include "merge-segments/merge-segments-main.h"

#include <fstream>

#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/util/file-functions.h>

#include "copy-number-caller/copy-number-caller-options.h"
#include "io/copy-number-caller-default-filenames.h"
#include "io/write-segments.h"
#include "likelihood/flag-copy-number-calls.h"
#include "mapq-utils.h"
#include "merge-segments/merge-segments.h"
#include "seg-to-vcf/seg-to-vcf.h"
#include "segmentation/read-segments.h"
#include "segmentation/segment-type.h"
#include "segmentation/segments-header.h"

namespace xoos::cnc {
using segmentation::SegmentType;

void MergeSegmentsMain(const CopyNumberCallerOptions& options) {
  const auto vcf_out_fname = options.output_dir / kDefaultGermlineCNCallsetVcfOutput;
  const auto merged_segments_out_fname = options.output_dir / kDefaultGermlineCNCallsetSegOutput;
  file::CheckFilePermissionsAndOutputPathExistence({}, {merged_segments_out_fname, vcf_out_fname});
  std::vector<GenomicSegment> segments =
      segmentation::ReadSegments(options.segments_fname.value(), SegmentType::kGermlineLikelihood);
  const segmentation::SegmentsHeader header = segmentation::ReadHeaderFromSegments(options.segments_fname.value());
  Observations mapqs;
  if (options.merge_segments_options.recalculate_per_segment_data ||
      options.merge_segments_options.use_mapqs_observations) {
    if (!options.mapping_qualities_fname.has_value()) {
      Logging::Error("Must provide --mapqs option!");
      throw std::runtime_error("missing option");
    }
    std::ifstream mapqs_ifs(options.mapping_qualities_fname.value());
    mapqs = ReadObservations(mapqs_ifs);
    if (options.merge_segments_options.use_mapqs_observations) {
      AssignAvgMeanMapqPerSegment(segments, mapqs);
    }
  }
  auto merged_segments = MergeAdjacentEqualCopyNumberSegments(segments);
  merged_segments = MergeLowMapqSegments(merged_segments, options.merge_segments_options.min_mapq_threshold);
  if (options.merge_segments_options.recalculate_per_segment_data) {
    if (!options.logrs_fname.has_value()) {
      Logging::Error("Must provide --logrs with the --recalculate-per-segment-data option!");
      throw std::runtime_error("missing option");
    }
    std::ifstream logrs_ifs(options.logrs_fname.value());
    Observations logrs = ReadObservations(logrs_ifs);
    AssignAvgMeanMapqPerSegment(merged_segments, mapqs);
    segmentation::PopulateGenomicSegmentOptionalFields(merged_segments, logrs);
    FlagCallsByExpectedTotalCopyNumber(merged_segments, options.sample_metadata_options.sex.value());
  }
  FlagCallsByMeanMapq(merged_segments, options.merge_segments_options.min_mapq_threshold);
  for (auto& seg : merged_segments) {
    if (!seg.flags.has_value()) {
      seg.flags.emplace(std::vector<std::string>{});
    }
  }
  WriteSegments(merged_segments_out_fname,
                merged_segments,
                segmentation::SegmentType::kGermlineLikelihood,
                options.command_line_info);
  WriteSegmentsToVcf(vcf_out_fname,
                     merged_segments,
                     {.seq_lengths = {},
                      .seq_order = {},
                      .sample_id = options.sample_metadata_options.sample_id,
                      .reference_file = header.reference_file,
                      .sex = options.sample_metadata_options.sex.value(),
                      .command_line_info = options.command_line_info});
}

}  // namespace xoos::cnc
