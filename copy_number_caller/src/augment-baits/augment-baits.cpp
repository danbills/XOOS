#include "augment-baits/augment-baits.h"

#include <algorithm>
#include <utility>

#include <fmt/core.h>

#include <xoos/log/logging.h>
#include <xoos/util/file-functions.h>

#include "baits.h"
#include "copy-number-caller/common/check-chr-y-par.h"
#include "copy-number-caller/copy-number-caller-options.h"
#include "intervals/interval-tools.h"
#include "intervals/intervals.h"
#include "io/copy-number-caller-default-filenames.h"
#include "segmentation/read-segments.h"
#include "segmentation/segment-type.h"
#include "sex.h"

namespace xoos::cnc {

void AugmentBaitsMain(const CopyNumberCallerOptions& options) {
  BaitRecords baits;
  const auto augmented_baits_out = options.output_dir / kDefaultAugmentedBaitsOutput;
  file::CheckFilePermissionsAndOutputPathExistence({}, {augmented_baits_out});
  if (options.augment_baits_options.generate_whole_genome_baits) {
    baits = GenerateAndAugmentWholeGenomeIntervals(options);
  } else {
    if (!options.augmented_baits_fname.has_value()) {
      Logging::Error("Baits file must be provided if not generating whole-genome intervals");
      throw std::runtime_error("Baits file must be provided if not generating whole-genome intervals");
    }
    baits = GenerateAndAugmentOnAndOffTargetBaits(options);
  }
  const vec<segmentation::GenomicSegment> seed_segments =
      segmentation::ReadSegments(options.seed_segments_fname, segmentation::SegmentType::kSeed);
  baits = RemovePARIfChrYPARUnmasked(options.normal_bam_fname, baits, seed_segments);
  std::ofstream ofs(augmented_baits_out);
  baits.Write(ofs, options.command_line_info);
}

BaitRecords GenerateAndAugmentWholeGenomeIntervals(const CopyNumberCallerOptions& options) {
  SeqLenMap seq_lengths(ParseFai(options.reference_genome_fai_fname));

  // Check if mappability file exists and is not empty
  if (!fs::exists(options.mappability_bigwig_fname)) {
    Logging::Error("Mappability file does not exist");
    throw std::runtime_error("Mappability file does not exist: " + options.mappability_bigwig_fname.string());
  }

  if (fs::file_size(options.mappability_bigwig_fname) == 0) {
    Logging::Error("Mappability file is empty");
    throw std::runtime_error("Empty mappability file: " + options.mappability_bigwig_fname.string());
  }

  BigWig mapp_bigwig(options.mappability_bigwig_fname);
  ContigToIntervals intervals;
  f64 default_min_mappability;
  Logging::Info("Generating {}-bp intervals across whole genome",
                options.augment_baits_options.whole_genome_interval_size);
  intervals = IntervalsFromSeqLenMap(seq_lengths, options.augment_baits_options.include_alt_contigs);
  ThrowIfEmpty(intervals);
  Logging::Info("{} intervals", GetNContigToIntervals(intervals));
  Logging::Info("Removing blocklisted regions from intervals");
  if (options.blocklist_bed_fname.has_value()) {
    const ContigToIntervals blocklist_intervals = BedToMap(options.blocklist_bed_fname.value());
    intervals = SubtractRegionsFromIntervals(intervals, blocklist_intervals);
    Logging::Info("{} intervals after removing blocklisted regions", GetNContigToIntervals(intervals));
  }
  default_min_mappability = options.augment_baits_options.whole_genome_min_mappability;
  Logging::Info("{} intervals", GetNContigToIntervals(intervals));
  Logging::Info("Cutting mappability=0 regions from intervals");
  intervals = Remove0MappableRegions(intervals, mapp_bigwig);
  ThrowIfEmpty(intervals, "No usable intervals. Check if mappability values overlap with reference genome");
  Logging::Info("{} intervals after removing 0-mappability regions", GetNContigToIntervals(intervals));
  Logging::Info("Splitting intervals into sub-intervals of length {}",
                options.augment_baits_options.whole_genome_interval_size);
  intervals = SplitIntervalsByWindowSizedLastWindowSizeOfRemainder(
      intervals, options.augment_baits_options.whole_genome_interval_size);
  ThrowIfEmpty(intervals);
  Logging::Info("{} intervals", GetNContigToIntervals(intervals));
  Logging::Info("Removing low mappability autosomal on-target intervals");
  intervals = RemoveLowMappableRegions(intervals, mapp_bigwig, default_min_mappability, SexChromHandling::kNoSexChrom);
  ThrowIfEmpty(intervals);
  Logging::Info("Removing low mappability allosomal on-target intervals");
  intervals = RemoveLowMappableRegions(intervals,
                                       mapp_bigwig,
                                       options.augment_baits_options.sex_chromosome_min_mappability,
                                       SexChromHandling::kOnlySexChrom);
  ThrowIfEmpty(intervals);
  Logging::Info("{} intervals", GetNContigToIntervals(intervals));
  // now put intervals into a Baits class
  BaitRecords baits;
  baits.SetRegions(ToRegionStringVec(intervals));
  baits.SetOnTargetStatus(std::vector<bool>(GetNContigToIntervals(intervals), true));
  // TODO: any modification of the intervals should be pushed earlier
  // the purpose of the Baits class should purely be for annotation
  intervals = RemoveRegionsByLength(intervals, options.augment_baits_options.whole_genome_interval_size);
  baits.SetRegions(ToRegionStringVec(intervals));  // Update baits with modified intervals
  baits.SetOnTargetStatus(std::vector<bool>(GetNContigToIntervals(intervals), true));  // Update on-target status
  baits.SetSeqLengths(seq_lengths);  // annotate with original contig lengths
  Logging::Info("Sorting regions");
  baits.SortByRegion(options.reference_genome_fai_fname);  // sort
  io::FastaReader fa(options.reference_genome_fname);

  Logging::Info("Annotating baits with GC bias");
  baits.AnnotateGC(fa);
  Logging::Info("Removing baits with extreme GC content");
  baits.FilterExtremeGC(options.augment_baits_options.min_gc_content, options.augment_baits_options.max_gc_content);

  Logging::Info("Annotating baits with mappability");
  baits.AnnotateMappability(mapp_bigwig);
  baits.SetReferenceFile(options.reference_genome_fname);
  return baits;
}

BaitRecords GenerateAndAugmentOnAndOffTargetBaits(const CopyNumberCallerOptions& options) {
  SeqLenMap seq_lengths(ParseFai(options.reference_genome_fai_fname));
  BigWig mapp_bigwig(options.mappability_bigwig_fname);
  ContigToIntervals on_targets;
  // if we're generating from whole genome, we also have to disable off-targets
  f64 default_min_mappability;
  Logging::Info("Generating on-targets from bed file");
  on_targets = BedToMap(options.augmented_baits_fname.value());
  ThrowIfEmpty(on_targets);
  Logging::Info("{} intervals", GetNContigToIntervals(on_targets));
  ContigToIntervals blocklist_intervals;
  if (options.blocklist_bed_fname.has_value()) {
    Logging::Info("Removing blocklisted regions from on-target intervals");
    blocklist_intervals = BedToMap(options.blocklist_bed_fname.value());
    on_targets = SubtractRegionsFromIntervals(on_targets, blocklist_intervals);
    Logging::Info("{} on-target intervals after removing blocklisted regions", GetNContigToIntervals(on_targets));
  }
  default_min_mappability = options.augment_baits_options.on_target_min_mappability;
  Logging::Info("{} on-target intervals", GetNContigToIntervals(on_targets));
  Logging::Info("Cutting mappability=0 regions from on-targets");
  on_targets = Remove0MappableRegions(on_targets, mapp_bigwig);
  ThrowIfEmpty(on_targets);
  Logging::Info("{} on-target intervals after removing 0-mappability regions", GetNContigToIntervals(on_targets));
  Logging::Info("Splitting on-target intervals into sub-intervals of average length {}",
                options.augment_baits_options.on_target_interval_size);
  on_targets = SplitIntervalsByWindowSizeAdjustForEqualWindows(on_targets,
                                                               options.augment_baits_options.on_target_interval_size);
  ThrowIfEmpty(on_targets);
  Logging::Info("{} on-target intervals", GetNContigToIntervals(on_targets));
  Logging::Info("Removing low mappability autosomal on-target intervals");
  on_targets =
      RemoveLowMappableRegions(on_targets, mapp_bigwig, default_min_mappability, SexChromHandling::kNoSexChrom);
  ThrowIfEmpty(on_targets);
  Logging::Info("Removing low mappability allosomal on-target intervals");
  on_targets = RemoveLowMappableRegions(on_targets,
                                        mapp_bigwig,
                                        options.augment_baits_options.sex_chromosome_min_mappability,
                                        SexChromHandling::kOnlySexChrom);
  ThrowIfEmpty(on_targets);
  Logging::Info("{} on-target intervals", GetNContigToIntervals(on_targets));
  // TODO: any modification of the intervals should be pushed earlier
  // the purpose of the Baits class should purely be for annotation
  BaitRecords baits;
  if (!options.augment_baits_options.no_off_targets) {
    Logging::Info("Generating off-target intervals");
    const ContigToIntervals off_targets =
        GenerateOffTargets(on_targets,
                           seq_lengths,
                           mapp_bigwig,
                           blocklist_intervals,
                           options.augment_baits_options.off_target_min_width,
                           options.augment_baits_options.off_target_interval_size,
                           options.augment_baits_options.off_target_trim_length,
                           options.augment_baits_options.off_target_min_mappability,
                           options.augment_baits_options.sex_chromosome_min_mappability);
    // concat on- and off- target intervals. they will be sorted within the baits object
    std::vector<std::string> on_target_regions = ToRegionStringVec(on_targets);
    std::vector<std::string> off_target_regions = ToRegionStringVec(off_targets);
    size_t n_on_target = on_target_regions.size();
    size_t n_off_target = off_target_regions.size();
    std::vector<std::string> all_regions;
    all_regions.reserve(n_on_target + n_off_target);
    std::move(on_target_regions.begin(), on_target_regions.end(), std::back_inserter(all_regions));
    std::move(off_target_regions.begin(), off_target_regions.end(), std::back_inserter(all_regions));
    // create the on-target std::vector
    std::vector<bool> on_target_status;
    on_target_status.reserve(all_regions.size());
    for (size_t i = 0; i < n_on_target; ++i) {
      on_target_status.push_back(true);
    }
    for (size_t i = 0; i < n_off_target; ++i) {
      on_target_status.push_back(false);
    }
    baits.SetRegions(std::move(all_regions));
    baits.SetOnTargetStatus(std::move(on_target_status));
  } else {
    Logging::Info("Off target generation was disabled");
    baits.SetRegions(ToRegionStringVec(on_targets));
    baits.SetOnTargetStatus(std::vector<bool>(GetNContigToIntervals(on_targets), true));
  }
  baits.SetSeqLengths(seq_lengths);
  Logging::Info("Sorting regions");
  baits.SortByRegion(options.reference_genome_fai_fname);
  io::FastaReader fa(options.reference_genome_fname);
  Logging::Info("Annotating baits with GC bias");
  baits.AnnotateGC(fa);
  Logging::Info("Annotating baits with mappability");
  baits.AnnotateMappability(mapp_bigwig);
  baits.SetReferenceFile(options.reference_genome_fname);
  return baits;
}
}  // namespace xoos::cnc
