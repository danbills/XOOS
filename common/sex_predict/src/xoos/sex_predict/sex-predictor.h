#pragma once

#include <filesystem>
#include <optional>
#include <vector>

#include <xoos/io/alignment-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "sex.h"

namespace xoos::sex_predict {

namespace fs = std::filesystem;

constexpr f32 kMinRatioForSexDet = 25;
constexpr f32 kMinRatioForSexNA = 20;

struct RegionInfo {
  s32 tid;
  std::string name;
};

std::optional<RegionInfo> GetXChromosome(const io::SamHdrPtr& bam_header);
std::optional<RegionInfo> GetYChromosome(const io::SamHdrPtr& bam_header);
std::optional<RegionInfo> GetChromosomeForStrings(const io::SamHdrPtr& bam_header,
                                                  const std::vector<std::string>& chromosomes);
Sex PredictSex(const fs::path& bam_file_path);

}  // namespace xoos::sex_predict
