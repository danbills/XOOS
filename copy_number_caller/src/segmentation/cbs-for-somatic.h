#pragma once

#include <armadillo>
#include <vector>

#include <xoos/enum/enum-util.h>
#include <xoos/types/float.h>

#include "segmentation/genomic-segments.h"

namespace xoos::cnc::segmentation::somatic {

f64 TStat(f64 n_v1, f64 mean_v1, f64 var_v1, f64 n_v2, f64 mean_v2, f64 var_v2);
f64 TStatFromVecs(const arma::vec& v1, const arma::vec& v2);
Segment GetMaxTSegCumSum(const arma::vec& obvs, size_t max_n, size_t min_n);
bool CircularPermuteNullAndTestSignificance(const arma::vec& obvs,
                                            size_t n_permutations,
                                            f64 max_p,
                                            f64 t,
                                            size_t max_obs_size,
                                            const arma::uvec& boundaries,
                                            size_t min_obs_size);
bool BinaryPermuteNullAndTestSignificance(
    const arma::vec& obvs, size_t n_permutations, f64 max_p, const Segment& seg1, const Segment& seg2);
f64 SiegmundPNuFull(f64 x);
f64 SiegmundP(f64 t, f64 m, f64 k);
bool HybridTestForSignificance(
    const arma::vec& obvs, size_t n_permutations, f64 max_p, f64 t, const arma::uvec& boundaries, size_t min_obs_size);
std::vector<Segment> GetSegments(
    const arma::vec& obvs, f64 max_p, size_t n_permutations, const arma::uvec& boundaries, size_t min_obs_size);
std::vector<Segment> CircularBinarySegmentation(
    const arma::vec& obvs, f64 max_p, size_t n_permutations, size_t min_obs_size, bool no_single_obvs_subsegments);
std::vector<Segment> MergeSingleObvsEndSegment(const std::vector<Segment>& sub_segments, size_t min_obs_size);
}  // namespace xoos::cnc::segmentation::somatic
