#pragma once
#include <armadillo>
#include <vector>

#include <xoos/types/float.h>

#include "observations.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc::segmentation {
f64 GetTrimmedStandardDeviation(const arma::vec& arr);
/**
 * @brief Undo genomic segments by merging adjacent segments if their median observation difference is less than a
 * factor (max_trimmed_sd_factor) of the trimmed standard deviation (trimmed_sd)
 * @param segments input segments
 * @param obvs observations vector
 * @param trimmed_sd trimmed standard deviation of the observations (use GetTrimmedStandardDeviation to calculate)
 * @param max_trimmed_sd_factor maximum factor of trimmed standard deviation to consider for merging segments
 * @return merged segments
 */
/// Merge two Observations by concatenating all fields in order: [obvs1, obvs2].
Observations MergeObservations(const Observations& obvs1, const Observations& obvs2);

std::vector<GenomicSegment> UndoSegments(const std::vector<GenomicSegment>& segments,
                                         const Observations& obvs,
                                         f64 trimmed_sd,
                                         f64 max_trimmed_sd_factor);
std::vector<GenomicSegment> IterativeUndoGenomicSegments(const std::vector<GenomicSegment>& seed_segments,
                                                         const std::vector<GenomicSegment>& segments,
                                                         const Observations& obvs,
                                                         f64 start_max_trimmed_sd_factor,
                                                         f64 increment_max_trimmed_sd_factor,
                                                         size_t max_num_segments);
}  // namespace xoos::cnc::segmentation
