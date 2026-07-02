#include "check-chr-y-par.h"

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "utility/utility-functions.h"

namespace xoos::cnc {
/**
 * @brief Checks whether the chrY PAR region was masked during alignment by looking for any alignments in the chrY PAR
 * region defined in the seed segments file. If any alignments are found in the chrY PAR region, an error is thrown
 * since this could significantly affect copy number predictions in allosomes. If the seed segments do not specify any
 * segments in the chrY PAR region, a warning is logged since we cannot check whether the chrY PAR region was masked
 * during alignment.
 */
bool CheckChrYPARUnmasked(const fs::path& bam_file, const vec<segmentation::GenomicSegment>& seed_segments) {
  // find the chrY seed segment that is in PAR region
  auto result = std::find_if(seed_segments.begin(), seed_segments.end(), [](const segmentation::GenomicSegment& seg) {
    return ((seg.contig == "chrY" || seg.contig == "Y") && seg.in_pseudo_autosomal_region.has_value() &&
            seg.in_pseudo_autosomal_region.value());
  });
  // no need to throw error if there aren't any seed segments in chrY PAR region, as that likely means the user is not
  // using any seed segments in the chrY PAR region for calling (which is reasonable since those regions can be noisy)
  if (result == seed_segments.end()) {
    Logging::Warn(
        "Could not find any seed segments in chrY PAR region. Cannot check whether chrY PAR region was masked during "
        "alignment. This could affect copy number predictions in allosomes.");
  } else {
    // open the BAM file and check for reads that are in the chrY PAR region described by the seed segment file.
    io::SamFilePtr sam_file = io::SamOpen(bam_file.string(), "r");
    io::SamHdrPtr sam_hdr = io::SamHdrRead(sam_file.get());
    const fs::path bai_file = bam_file.string() + ".bai";
    io::HtsIdxPtr sam_idx = io::SamIndexLoad(sam_file.get(), bai_file);
    std::string region_string =
        result->contig + ":" + std::to_string(result->start + 1) + "-" + std::to_string(result->end);
    io::HtsItrPtr sam_itr_ptr = io::SamItrQuerySNoThrow(sam_idx.get(), sam_hdr.get(), region_string);
    // if null is returned by sam_itr_query, then it means that chrY or the PAR region was not in the reference FASTA
    // used by the aligner. Not issuing a warning because it is not deletrious behavior.
    if (sam_itr_ptr.get() == nullptr) {
      Logging::Info("BAM file has no alignments in the chrY PAR regions");
    } else {
      size_t n_alignments = 0;
      const io::Bam1Ptr bam1_ptr(bam_init1());
      while (sam_itr_next(sam_file.get(), sam_itr_ptr.get(), bam1_ptr.get()) > 0) {
        // only increment n_alignments if the read starts AFTER the start of the PAR region and ends BEFORE the end of
        // the PAR region, since there can be alignments that partially overlap the PAR region but are not fully
        // contained within it
        if (bam1_ptr->core.pos >= ToSigned(result->start) && bam_endpos(bam1_ptr.get()) <= ToSigned(result->end)) {
          ++n_alignments;
        }
      }
      if (n_alignments > 0) {
        Logging::Warn(
            "Found {} alignments in chrY PAR region. This suggests that the chrY PAR region may NOT have been masked "
            "during alignment. This could affect copy number predictions in allosomes.",
            n_alignments);
        return true;
      }
    }
  }
  return false;
}

BaitRecords RemovePARFromBaits(const BaitRecords& baits, const vec<segmentation::GenomicSegment>& seed_segments) {
  // Get the PAR regions from the seed segments
  vec<std::tuple<std::string, size_t, size_t>> par_regions;
  for (const auto& seg : seed_segments) {
    if (seg.in_pseudo_autosomal_region.has_value() && seg.in_pseudo_autosomal_region.value()) {
      par_regions.emplace_back(seg.contig, seg.start, seg.end);
    }
  }
  if (par_regions.empty()) {
    Logging::Warn("No seed segments found in chrY PAR region. Cannot remove PAR regions from augmented baits.");
    return baits;
  }
  // Remove baits that overlap with the PAR regions
  std::vector<std::string> filtered_regions;
  vec<f64> filtered_gc_bias;
  vec<f64> filtered_mappability;
  std::vector<bool> filtered_on_target;
  const auto& regions = baits.GetRegions();
  const auto& gc_biases = baits.GetGCBias();
  const auto& mappabilities = baits.GetMappability();
  const auto& on_targets = baits.GetOnTargetStatus();
  for (size_t i = 0; i < baits.GetRegions().size(); ++i) {
    // parse the region string to get the contig, start, and end
    auto [contig, start, end] = ParseRegionString(regions[i]);
    bool in_par_region = std::any_of(par_regions.begin(),
                                     par_regions.end(),
                                     [&contig, start, end](const std::tuple<std::string, size_t, size_t>& par_region) {
                                       const auto& [par_contig, par_start, par_end] = par_region;
                                       return (contig == par_contig && (start <= par_end && end >= par_start));
                                     });
    if (!in_par_region) {
      filtered_regions.push_back(regions[i]);
      filtered_gc_bias.push_back(gc_biases[i]);
      filtered_mappability.push_back(mappabilities[i]);
      filtered_on_target.push_back(on_targets[i]);
    }
  }
  BaitRecords filtered_baits;
  filtered_baits.SetRegions(std::move(filtered_regions));
  // the following two are not using std::move because clangd-tidy produces the "performance-move-const-arg" warning
  // probably related to armadillo?
  filtered_baits.SetGCBias(filtered_gc_bias);
  filtered_baits.SetMappability(filtered_mappability);
  filtered_baits.SetOnTargetStatus(std::move(filtered_on_target));
  auto seq_lengths = baits.GetSeqLengths();
  filtered_baits.SetSeqLengths(seq_lengths);
  filtered_baits.SetReferenceFile(baits.GetReferenceFile());
  Logging::Info(
      "Removed baits that overlap with PAR regions. Original number of baits: {}, number of baits after filtering: {}.",
      baits.GetRegions().size(),
      filtered_baits.GetRegions().size());
  return filtered_baits;
}

BaitRecords RemovePARIfChrYPARUnmasked(const fs::path& bam_file,
                                       const BaitRecords& baits,
                                       const vec<segmentation::GenomicSegment>& seed_segments) {
  if (CheckChrYPARUnmasked(bam_file, seed_segments)) {
    return RemovePARFromBaits(baits, seed_segments);
  }
  Logging::Info("chrY PAR region appears to be masked during alignment. PAR regions will be included in analysis");
  return baits;
}
}  // namespace xoos::cnc
