#pragma once
#include <memory>
#include <optional>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>

#include "augment-baits/augment-baits-options.h"
#include "calculate-coverage/calculate-coverage-options.h"
#include "copy-number-caller/copy-number-caller-modes.h"
#include "denoise/denoise-options.h"
#include "gc-correct/gc-correct-options.h"
#include "likelihood/likelihood-options.h"
#include "merge-segments/merge-segments-options.h"
#include "misc/sample-metadata-options.h"
#include "misc/vcf-parsing-options.h"
#include "purity-ploidy-search/purity-ploidy-search-options.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc {

const size_t kCopyNumberCallerDefaultTumorMinCoverage = 0;
const size_t kCopyNumberCallerDefaultSomaticNormalMinCoverage = 30;
const size_t kCopyNumberCallerDefaultGermlineNormalMinCoverage = 0;

struct CopyNumberCallerOptions {
  // Algorithm Submodule Options
  AugmentBaitsOptions augment_baits_options{};
  CalculateCoverageOptions calculate_coverage_options{};
  GCCorrectOptions gc_correct_options{};
  DenoiseOptions denoise_options{};
  SegmentationOptions segmentation_options{};
  PurityPloidySearchOptions purity_ploidy_search_options{};
  LikelihoodOptions likelihood_options{};
  MergeSegmentsOptions merge_segments_options{};
  BafFilterOptions baf_filter_options{};
  SampleMetadataOptions sample_metadata_options{};

  size_t threads{1};
  CopyNumberCallerModes mode{CopyNumberCallerModes::kUnknown};

  // Somatic Coverage Thresholds
  size_t tumor_min_coverage{};
  size_t normal_min_coverage{};

  // Somatic Input File
  fs::path tumor_bam_fname{};
  std::optional<fs::path> tumor_coverage_fname{};

  // Input Files
  fs::path normal_bam_fname{};
  fs::path reference_genome_fname{};
  fs::path reference_genome_fai_fname{};
  fs::path mappability_bigwig_fname{};
  fs::path seed_segments_fname{};
  std::optional<fs::path> blocklist_bed_fname{};
  std::optional<fs::path> augmented_baits_fname{};
  std::optional<fs::path> normal_coverage_fname{};
  std::optional<fs::path> vcf_fname{};
  std::optional<fs::path> ad_fname{};
  std::optional<fs::path> segments_fname{};
  std::optional<fs::path> logrs_fname{};
  std::optional<fs::path> bafs_fname{};
  std::optional<fs::path> mapping_qualities_fname{};

  // Panel of Normals
  fs::path panel_of_normals_lists{};
  fs::path panel_of_normals_hdf5_fname{};

  // Output Directory
  fs::path output_dir{};

  // Output-control flags (hidden advanced options)
  bool save_taskflow_graph{false};
  bool save_log_ratio_segments{false};
  bool save_baf_segments{false};

  // Structured metadata for all output formats (TSV and VCF)
  io::CommandLineInfo command_line_info;
};

using CopyNumberCallerOptionsPtr = std::shared_ptr<CopyNumberCallerOptions>;

}  // namespace xoos::cnc
