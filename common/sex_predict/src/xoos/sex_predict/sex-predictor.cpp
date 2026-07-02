#include "sex-predictor.h"

#include <algorithm>
#include <filesystem>
#include <vector>

#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/io/alignment-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>

namespace xoos::sex_predict {

// get the most likely X chromosome string from the bam header
std::optional<RegionInfo> GetChromosomeForStrings(const io::SamHdrPtr& bam_header,
                                                  const std::vector<std::string>& chromosomes) {
  if (chromosomes.empty()) {
    return std::nullopt;
  }
  // iterate through all regions in the header and find the X chromosome
  for (s32 i = 0; i < bam_header->n_targets; ++i) {
    // convert bam_header->target_name[i] to a lower-case string
    std::string target_name = bam_header->target_name[i];
    std::ranges::transform(target_name, target_name.begin(), ::tolower);
    if (std::ranges::find(chromosomes, target_name) != chromosomes.end()) {
      return RegionInfo{i, bam_header->target_name[i]};
    }
  }
  return std::nullopt;
}

std::optional<RegionInfo> GetXChromosome(const io::SamHdrPtr& bam_header) {
  return GetChromosomeForStrings(bam_header, std::vector<std::string>{"chrx", "x"});
}

std::optional<RegionInfo> GetYChromosome(const io::SamHdrPtr& bam_header) {
  return GetChromosomeForStrings(bam_header, std::vector<std::string>{"chry", "y"});
}

// Get all the reads from the bam file index and calculate the xy ration between X and Y.
// If the ratio X / Y > 25 then return female
// If the ratio X / Y > 20 then return unknown
// Else return male.
Sex PredictSex(const fs::path& bam_file_path) {
  const auto [bam_file, header, idx] = io::OpenAlignmentReader(bam_file_path);
  const auto alignment_reference_regions = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>(bam_init1(), bam_destroy1);

  // get the number of reads for X and Y

  // count X reads
  const auto region_x = GetXChromosome(header);
  u64 count_x = 0;
  if (region_x.has_value()) {
    const auto iter_x = io::HtsItrMultiPtr(sam_itr_querys(idx.get(), header.get(), region_x.value().name.c_str()));
    if (!iter_x) {
      throw error::Error("Could not parse region '{}'", region_x.value().name);
    }
    u64 unmapped_x = 0;
    // use hts_idx_get_stat to get the number of reads in the region
    hts_idx_get_stat(idx.get(), region_x.value().tid, &count_x, &unmapped_x);
  }

  // count Y reads
  const auto region_y = GetYChromosome(header);
  u64 count_y = 0;
  if (region_y.has_value()) {
    const auto iter_y = io::HtsItrMultiPtr(sam_itr_querys(idx.get(), header.get(), region_y.value().name.c_str()));
    if (!iter_y) {
      throw error::Error("Could not parse region '{}'", region_y.value().name.c_str());
    }

    u64 unmapped_y = 0;
    hts_idx_get_stat(idx.get(), region_y.value().tid, &count_y, &unmapped_y);
  }

  // and calculate the ratio between X and Y
  const f64 xy_ratio = static_cast<f64>(count_x) / static_cast<f64>(count_y);
  if (std::isnan(xy_ratio)) {
    return Sex::kUnknown;
  }
  if (xy_ratio > kMinRatioForSexDet) {
    return Sex::kFemale;
  }
  if (xy_ratio > kMinRatioForSexNA) {
    return Sex::kUnknown;
  }
  return Sex::kMale;
}

}  // namespace xoos::sex_predict
