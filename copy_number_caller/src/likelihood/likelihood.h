#pragma once
#include <armadillo>
#include <cmath>

#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "copy-number-caller/copy-number-caller-options.h"
#include "likelihood/likelihood-mode.h"
#include "likelihood/likelihood-model.h"
#include "likelihood/likelihood-options.h"
#include "observations.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/interval-trees.h"

namespace xoos::cnc {
using segmentation::GenomicSegment;
using segmentation::IntervalTrees;

const f64 kNegativeInf = -1 * arma::datum::inf;

struct SegmentAllSummarizedBAFLikelihoodsResults {
  arma::vec baf_likelihoods;
};

void CalculateLikelihoodsMain(const CopyNumberCallerOptions& options);
std::vector<GenomicSegment>& CalculateLikelihoods(std::vector<GenomicSegment>& segments,
                                                  const Observations& logrs,
                                                  const std::optional<Observations>& ref_ads,
                                                  const std::optional<Observations>& alt_ads,
                                                  const std::optional<Observations>& mapqs,
                                                  LikelihoodModel ll_model,
                                                  const LikelihoodOptions& options,
                                                  const SampleMetadataOptions& sample_metadata_options);
std::vector<GenomicSegment>& AllSegmentsLogRLikelihoods(std::vector<GenomicSegment>& segments,
                                                        const Observations& logrs,
                                                        f64 purity,
                                                        f64 ploidy,
                                                        LikelihoodMode mode,
                                                        Sex sex);
std::vector<GenomicSegment>& AllSegmentsSerialLogRBAFLikelihoods(std::vector<GenomicSegment>& segments,
                                                                 const Observations& logrs,
                                                                 const Observations& ref_ads,
                                                                 const Observations& alt_ads,
                                                                 f64 purity,
                                                                 f64 ploidy,
                                                                 LikelihoodMode mode,
                                                                 Sex sex);
f64 SegmentSummarizedLogRLikelihood(
    const arma::vec& logrs, f64 sd, f64 purity, f64 tumor_ploidy, f64 total_copy_number, bool expect_haploid);
f64 SegmentSummarizedBAFLikelihood(const arma::vec& ref_ads,
                                   const arma::vec& alt_ads,
                                   f64 sd,
                                   f64 purity,
                                   f64 total_copy_number,
                                   f64 multiplicity,
                                   bool expect_haploid);
f64 GetMeanSegmentLogRSD(const std::vector<GenomicSegment>& segments,
                         const Observations& logrs,
                         IntervalTrees& logr_trees);
f64 GetMeanSegmentMBAFSD(const std::vector<GenomicSegment>& segments,
                         const Observations& ref_ads,
                         const Observations& alt_ads,
                         IntervalTrees& variant_trees);
arma::vec& SegmentAllSummarizedLogRLikelihoods(
    const arma::vec& logrs, f64 sd, f64 purity, f64 tumor_ploidy, bool expect_haploid, arma::vec& ret);
SegmentAllSummarizedBAFLikelihoodsResults SegmentAllSummarizedBAFLikelihoods(
    f64 mbaf, f64 sd, f64 purity, f64 total_copy_number, bool expect_haploid);
GenomicSegment& SegmentSerialSummarizedLogRBAFLikelihood(GenomicSegment& seg,
                                                         const arma::vec& logrs,
                                                         f64 purity,
                                                         f64 ploidy,
                                                         f64 mean_logr_sd,
                                                         f64 mean_baf_sd,
                                                         bool prior_2n,
                                                         bool expect_haploid,
                                                         arma::vec& logr_likelihoods,
                                                         arma::vec& baf_likelihoods);
bool SkipBasedOnModeAndSex(const GenomicSegment& seg, LikelihoodMode mode, Sex sex);
bool ExpectHaploid(const GenomicSegment& seg, LikelihoodMode mode, Sex sex);
}  // namespace xoos::cnc
