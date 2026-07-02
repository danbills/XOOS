#pragma once

#include <armadillo>
#include <vector>

#include <xoos/types/float.h>

#include "genomic-segments.h"
#include "segmentation/max-t.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc::segmentation {

bool CircularPermuteNullAndTestSignificance(const arma::vec& obvs,
                                            size_t n_permutations,
                                            f64 max_p,
                                            Segment max_t_seg,
                                            size_t max_obs_size,
                                            const arma::uvec& boundaries,
                                            size_t min_obs_per_segment,
                                            CbsMaxTMethod max_t_method = CbsMaxTMethod::kBruteForce);
bool BinaryPermuteNullAndTestSignificance(const arma::vec& obvs,
                                          size_t n_permutations,
                                          f64 max_p,
                                          const Segment& seg1,
                                          const Segment& seg2,
                                          f64 min_t_for_automatic_reject_null);
f64 SiegmundPNuFull(f64 x);
f64 SiegmundP(f64 t, f64 m, f64 k);
bool HybridTestForSignificance(const arma::vec& obvs,
                               size_t n_permutations,
                               f64 max_p,
                               Segment max_t_seg,
                               const arma::uvec& boundaries,
                               size_t min_obs_per_segment,
                               CbsMaxTMethod max_t_method = CbsMaxTMethod::kBruteForce,
                               f64 min_t_for_automatic_reject_null = kSegmentationDefaultMinTForAutomaticSegmentation);
std::vector<Segment> GetSegments(
    const arma::vec& obvs,
    f64 max_p,
    size_t n_permutations,
    const arma::uvec& boundaries,
    size_t min_obs_per_segment,
    CbsMaxTMethod max_t_method = CbsMaxTMethod::kBruteForce,
    f64 min_t_for_automatic_reject_null = kSegmentationDefaultMinTForAutomaticSegmentation);
std::vector<Segment> CircularBinarySegmentation(const arma::vec& obvs, const SegmentationOptions& options);
std::vector<Segment> MergeSingleObvsEndSegment(const std::vector<Segment>& sub_segments, size_t min_obs_per_segment);
arma::vec TruncateOutliers(const arma::vec& obvs);

}  // namespace xoos::cnc::segmentation
