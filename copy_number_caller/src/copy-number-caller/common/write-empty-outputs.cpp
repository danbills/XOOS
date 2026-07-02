#include "copy-number-caller/common/write-empty-outputs.h"

#include <fstream>

#include <csv.hpp>

#include <xoos/io/metadata-util.h>

#include "coverage.h"
#include "io/column-names.h"
#include "io/write-bigwig.h"
#include "io/write-igv-xml.h"
#include "io/write-metrics.h"
#include "io/write-segments.h"
#include "mapq-utils.h"
#include "observations.h"
#include "purity-ploidy-search/purity-ploidy-search.h"
#include "seg-to-vcf/seg-to-vcf.h"

namespace xoos::cnc {

void WriteGermlineEmptyOutputs(const GermlineNormalWGSOutputPaths& paths,
                               const CoverageCheckResult& coverage_check,
                               const BaitRecords& baits,
                               const CopyNumberCallerOptions& options) {
  // SEG file — header only, no data rows
  WriteSegments(paths.likelihood_out, {}, segmentation::SegmentType::kGermlineLikelihood, options.command_line_info);

  // VCF file — header only, no records
  segmentation::WriteSegmentsToVcf(paths.vcf_out,
                                   {},
                                   {.seq_lengths = baits.GetSeqLengths(),
                                    .seq_order = GetContigOrder(options.reference_genome_fai_fname),
                                    .sample_id = options.sample_metadata_options.sample_id,
                                    .reference_file = baits.GetReferenceFile(),
                                    .sex = Sex::kUnknown,
                                    .mode = options.mode,
                                    .purity = 0.99,
                                    .ploidy = 2.0,
                                    .command_line_info = options.command_line_info});

  // Metrics file
  WriteMetricsFile(paths.metrics_out,
                   baits,
                   0,
                   Sex::kUnknown,
                   {},
                   std::nullopt,
                   std::nullopt,
                   options.mode,
                   options.command_line_info,
                   coverage_check);

  // Corrected coverage — header only
  {
    std::ofstream ofs(paths.normal_corrected_coverage_out);
    CoverageRecords{}.Write(ofs, options.command_line_info);
  }

  // Log ratios — header only
  {
    std::ofstream ofs(paths.logrs_out);
    Observations{}.Write(ofs, kColumnLogRatio, false, options.command_line_info);
  }

  // Log ratios bigwig
  WriteEmptyBigWig(paths.logrs_bw_out, options.reference_genome_fai_fname);

  // Mapping quality — header only
  {
    std::ofstream ofs(paths.mapping_qualities_out);
    Observations{}.Write(ofs, kMeanMapqColumnName, false, options.command_line_info);
  }

  // Log ratio segments (optional)
  if (paths.logr_segments_out.has_value()) {
    WriteSegments(paths.logr_segments_out.value(), {}, segmentation::SegmentType::kLogROnly, options.command_line_info);
  }

  // BAF files (optional)
  if (paths.baf_out.has_value()) {
    std::ofstream ofs(paths.baf_out.value());
    Observations{}.Write(ofs, kColumnBAF, false, options.command_line_info);
  }
  if (paths.baf_bw_out.has_value()) {
    WriteEmptyBigWig(paths.baf_bw_out.value(), options.reference_genome_fai_fname);
  }

  // IGV XML
  WriteIGVXML(paths.igv_xml_out, std::make_optional(paths.logrs_bw_out), paths.baf_bw_out);
}

void WriteSomaticEmptyOutputs(const SomaticTumorNormalWGSOutputPaths& paths,
                              const std::optional<CoverageCheckResult>& coverage_check,
                              const BaitRecords& baits,
                              const CopyNumberCallerOptions& options,
                              const std::optional<VcfCheckResult>& vcf_check) {
  // SEG file — header only
  WriteSegments(
      paths.likelihood_out, {}, segmentation::SegmentType::kSomaticWithBafLikelihood, options.command_line_info);

  // VCF file — header only
  segmentation::WriteSegmentsToVcf(paths.vcf_out,
                                   {},
                                   {.seq_lengths = baits.GetSeqLengths(),
                                    .seq_order = GetContigOrder(options.reference_genome_fai_fname),
                                    .sample_id = options.sample_metadata_options.sample_id,
                                    .reference_file = baits.GetReferenceFile(),
                                    .sex = Sex::kUnknown,
                                    .mode = options.mode,
                                    .purity = 0.0,
                                    .ploidy = 2.0,
                                    .vcf_purity_source = VcfPuritySource::kNone,
                                    .command_line_info = options.command_line_info});

  // Metrics file
  WriteMetricsFile(paths.metrics_out,
                   baits,
                   0,
                   Sex::kUnknown,
                   {},
                   std::nullopt,
                   std::nullopt,
                   options.mode,
                   options.command_line_info,
                   coverage_check,
                   vcf_check);

  // Corrected coverage — header only for both tumor and normal
  {
    std::ofstream ofs(paths.tumor_corrected_coverage_out);
    CoverageRecords{}.Write(ofs, options.command_line_info);
  }
  {
    std::ofstream ofs(paths.normal_corrected_coverage_out);
    CoverageRecords{}.Write(ofs, options.command_line_info);
  }

  // Log ratios — header only
  {
    std::ofstream ofs(paths.logrs_out);
    Observations{}.Write(ofs, kColumnLogRatio, false, options.command_line_info);
  }

  // Log ratios bigwig
  WriteEmptyBigWig(paths.logrs_bw_out, options.reference_genome_fai_fname);

  // Log ratio segments — header only
  WriteSegments(paths.logr_segments_out, {}, segmentation::SegmentType::kLogROnly, options.command_line_info);

  // BAF segments (optional)
  if (paths.baf_segments_out.has_value()) {
    WriteSegments(paths.baf_segments_out.value(), {}, segmentation::SegmentType::kBaf, options.command_line_info);
  }

  // Mapping quality — header only
  {
    std::ofstream ofs(paths.mapping_qualities_out);
    Observations{}.Write(ofs, kMeanMapqColumnName, false, options.command_line_info);
  }

  // Purity/ploidy grid — only when purity/ploidy are not user-supplied
  if (!options.sample_metadata_options.purity.has_value() || !options.sample_metadata_options.ploidy.has_value()) {
    std::ofstream ofs(paths.purity_ploidy_grid_out);
    if (!options.command_line_info.version.empty()) {
      io::WriteTsvMetadata(ofs, options.command_line_info);
    }
    ofs << kPurityPloidyGridOutPurity << "\t" << kPurityPloidyGridOutPloidy << "\t"
        << kPurityPloidyGridOutTotalJointLikelihood << "\t" << kPurityPloidyGridOutTotalLogRatioLikelihood << "\t"
        << kPurityPloidyGridOutTotalBAlleleFrequencyLikelihood << '\n';
  }

  // BAF files — header only
  {
    std::ofstream ofs(paths.baf_out);
    Observations{}.Write(ofs, kColumnBAF, false, options.command_line_info);
  }
  WriteEmptyBigWig(paths.baf_bw_out, options.reference_genome_fai_fname);

  // IGV XML
  WriteIGVXML(paths.igv_xml_out, std::make_optional(paths.logrs_bw_out), std::make_optional(paths.baf_bw_out));
}

}  // namespace xoos::cnc
