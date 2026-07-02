#include "likelihood/flag-copy-number-calls.h"

#include <xoos/log/logging.h>
#include <xoos/stats/copy-number-stats.h>

#include "likelihood/likelihood-flags.h"
#include "likelihood/likelihood.h"
#include "segmentation/interval-trees.h"
#include "sex.h"

namespace xoos::cnc {
using segmentation::IntervalTrees;

void FlagCallsByMeanMapq(std::vector<GenomicSegment>& segments, f64 mapq_cutoff_for_calls) {
  for (auto& segment : segments) {
    // skipped segments or segments with no MAPQ observations
    if (!segment.avg_mean_mapq.has_value()) {
      continue;
    }
    if ((!segment.total_copy_number.has_value() ||
         (segment.total_copy_number.has_value() && segment.total_copy_number.value() != 0)) &&
        segment.avg_mean_mapq.value() < mapq_cutoff_for_calls) {
      if (!segment.flags.has_value()) {
        segment.flags.emplace(std::vector<std::string>{});
      }
      segment.flags.value().emplace_back(kLikelihoodFlagLowAvgMeanMapq);
    }
  }
}

const f64 kLikelihoodMinTCNForHetDel = 0.7;
const f64 kLikelihoodMaxTCNForHetDel = 1.3;
const f64 kLikelihoodMaxTCNForHomDel = 0.3;

/**
 * @brief assigns the expected total copy number and then adds a flag if expected total copy number doesn't meet the
 * range for hom del or het del. This should only be used in germline mode. TODO break this up into two functions.
 * @param segments
 * @param logrs
 * @param sex
 */
void FlagCallsByExpectedTotalCopyNumber(std::vector<GenomicSegment>& segments, const Sex& sex) {
  for (auto& seg : segments) {
    bool expect_haploid = ExpectHaploid(seg, LikelihoodMode::kGermline, sex);
    f64 mean_logr = seg.mean_logr.value();
    f64 expected_total_copy_number =
        stats::ExpectedTotalCopyNumber(seg.purity.value(), seg.ploidy.value(), mean_logr, expect_haploid);
    seg.expected_total_copy_number = expected_total_copy_number;
    bool het_del_not_in_range =
        *seg.total_copy_number == 1 && (expected_total_copy_number < kLikelihoodMinTCNForHetDel ||
                                        expected_total_copy_number > kLikelihoodMaxTCNForHetDel);
    bool hom_del_not_in_range = *seg.total_copy_number == 0 && expected_total_copy_number > kLikelihoodMaxTCNForHomDel;
    auto flag_to_assign = kLikelihoodFlagNonIntegerTcn;
    // remove the flag if it's already there
    if (seg.flags.has_value()) {
      auto& flags = seg.flags.value();
      if (std::find(flags.begin(), flags.end(), flag_to_assign) != flags.end()) {
        flags.erase(std::remove(flags.begin(), flags.end(), flag_to_assign), flags.end());
      }
    }
    if (het_del_not_in_range || hom_del_not_in_range) {
      if (!seg.flags.has_value()) {
        seg.flags.emplace(std::vector<std::string>{});
      }
      seg.flags.value().emplace_back(flag_to_assign);
    }
  }
}

void FlagCallsByLength(std::vector<GenomicSegment>& segments, size_t cnv_length_flag_min_size) {
  for (auto& seg : segments) {
    if (seg.end - seg.start < cnv_length_flag_min_size) {
      if (!seg.flags.has_value()) {
        seg.flags.emplace(std::vector<std::string>{});
      }
      seg.flags.value().emplace_back(kLikelihoodFlagCnvLength);
    }
  }
}

/**
 * @brief adds the kLikelihoodFlagChrYFemale flag to the segment "flag" members that are on chrY if the sample is female
 * @param segments
 */
void FlagChrYCallsIfFemale(std::vector<GenomicSegment>& segments) {
  for (auto& seg : segments) {
    if (IsInChromY(seg.contig) && seg.sex == Sex::kFemale) {
      if (!seg.flags.has_value()) {
        seg.flags.emplace(std::vector<std::string>{});
      }
      seg.flags.value().emplace_back(kLikelihoodFlagChrYFemale);
    }
  }
}

}  // namespace xoos::cnc
