#pragma once
#include <armadillo>
#include <tuple>
#include <vector>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "copy-number-caller/copy-number-caller-options.h"
#include "likelihood/total-copy-number-prior.h"
#include "observations.h"
#include "segmentation/genomic-segments.h"
#include "sex.h"

using std::tuple;

namespace xoos::cnc {
using segmentation::GenomicSegment;

struct PurityPloidySearchOut {
  f64 purity = 0;
  f64 ploidy = 0;
  arma::mat joint_ll_grid;
  arma::mat logr_ll_grid;
  arma::mat baf_ll_grid;
  TotalCopyNumberPrior total_copy_number_prior = TotalCopyNumberPrior::kDiploidBias;
};

const std::string kPurityPloidyGridOutPurity = "Purity";
const std::string kPurityPloidyGridOutPloidy = "Ploidy";
const std::string kPurityPloidyGridOutTotalJointLikelihood = "TotalJointLikelihood";
const std::string kPurityPloidyGridOutTotalLogRatioLikelihood = "TotalLogRatioLikelihood";
const std::string kPurityPloidyGridOutTotalBAlleleFrequencyLikelihood = "TotalBAlleleFrequencyLikelihood";
const std::string kPurityPloidyOutPurity = "Purity";
const std::string kPurityPloidyOutPloidy = "Ploidy";

void PurityPloidySearchMain(const CopyNumberCallerOptions& options);
PurityPloidySearchOut PurityPloidySearch(const std::vector<GenomicSegment>& segments,
                                         const Observations& logrs,
                                         const std::optional<Observations>& ref_obvs,
                                         const std::optional<Observations>& alt_obvs,
                                         f64 seg_max_logr,
                                         size_t seg_min_num_logrs,
                                         size_t seg_min_num_snps,
                                         const fs::path& grid_out_fname,
                                         Sex sex,
                                         const io::CommandLineInfo& command_line_info,
                                         bool use_diploid_bias_prior);
PurityPloidySearchOut PurityPloidySearchWithWGDQC(const std::vector<GenomicSegment>& segments,
                                                  const Observations& logrs,
                                                  const std::optional<Observations>& ref_obvs,
                                                  const std::optional<Observations>& alt_obvs,
                                                  const PurityPloidySearchOptions& options,
                                                  const fs::path& grid_out_fname,
                                                  Sex sex,
                                                  const io::CommandLineInfo& command_line_info);
arma::mat GetPurityPloidyGridLogROnly(const std::vector<GenomicSegment>& segments,
                                      const Observations& obvs,
                                      f64 seg_max_logr,
                                      size_t seg_min_num_logrs,
                                      Sex sex,
                                      bool use_diploid_bias_prior = true);
std::tuple<arma::mat, arma::mat, arma::mat> GetPurityPloidyGridSerialLogRBAF(
    const std::vector<GenomicSegment>& segments,
    const Observations& logrs,
    const Observations& ref_ads,
    const Observations& alt_ads,
    f64 seg_max_logr,
    size_t seg_min_num_logrs,
    size_t seg_min_num_snps,
    Sex sex,
    bool use_diploid_bias_prior = true);
tuple<f64, f64> EstimatePurityPloidyByTotalMaximum(const arma::mat& purity_ploidy_grid);
tuple<f64, f64> EstimatePurityPloidyByLocalMaxima(const arma::mat& purity_ploidy_grid,
                                                  const Observations& ref_ads,
                                                  const Observations& alt_ads);
void WritePurityPloidyGrid(const arma::mat& purity_ploidy_grid_joint_lls,
                           const std::optional<arma::mat>& purity_ploidy_grid_logr_lls,
                           const std::optional<arma::mat>& purity_ploidy_grid_baf_lls,
                           const fs::path& out_fname,
                           const io::CommandLineInfo& command_line_info,
                           std::optional<TotalCopyNumberPrior> prior);
void WritePurityPloidy(f64 purity,
                       f64 ploidy,
                       const fs::path& out_fname,
                       const io::CommandLineInfo& command_line_info,
                       std::optional<TotalCopyNumberPrior> prior);
bool ShouldSkipSegment(const GenomicSegment& seg,
                       f64 seg_max_logr,
                       size_t seg_min_num_logrs,
                       std::optional<size_t> seg_min_num_snps,
                       bool check_seg_min_snps = true);
bool HasLowIntervalDensity(const std::vector<size_t>& logr_idxs, const Observations& logrs, const GenomicSegment& seg);

}  // namespace xoos::cnc
