#pragma once
#include <optional>
#include <string_view>
#include <vector>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "copy-number-caller/copy-number-caller-modes.h"
#include "copy-number-caller/copy-number-caller-options.h"
#include "io/fai.h"
#include "likelihood/total-copy-number-prior.h"
#include "segmentation/genomic-segments.h"
#include "sex.h"
#include "vcf-purity-source.h"

namespace xoos::cnc::segmentation {

inline constexpr std::string_view kSampleWarningHighExtremeBafProportion = "##SAMPLE_WARNING=highExtremeBafProportion";

struct WriteSegmentsToVcfParams {
  SeqLenMap seq_lengths;
  std::unordered_map<std::string, size_t> seq_order;
  std::string sample_id;
  std::string reference_file;
  Sex sex;
  CopyNumberCallerModes mode{CopyNumberCallerModes::kGermlineNormalWGS};
  std::optional<f64> purity;
  std::optional<f64> ploidy;
  std::optional<VcfPuritySource> vcf_purity_source;
  std::optional<TotalCopyNumberPrior> total_copy_number_prior;
  // Structured metadata for VCF ##RocheCommandLine header
  io::CommandLineInfo command_line_info;
  bool has_high_extreme_baf_proportion = false;
};

enum class CopyNumberStatus {
  kLoss,
  kRef,
  kCnLoh,
  kGain,
  kGainLoh,
  kUnknown
};

void WriteSegmentsToVcfMain(const CopyNumberCallerOptions& options);

void WriteSegmentsToVcf(const fs::path& fname,
                        const std::vector<GenomicSegment>& segments,
                        const WriteSegmentsToVcfParams& params);

CopyNumberStatus GetCopyNumberStatus(const GenomicSegment& seg);
CopyNumberStatus GetSomaticCopyNumberStatus(const GenomicSegment& seg);
std::string ConstructVcfFilterFieldFromGenomicSegment(const GenomicSegment& seg);
}  // namespace xoos::cnc::segmentation
