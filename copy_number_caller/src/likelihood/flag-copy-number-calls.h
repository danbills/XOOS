#pragma once
#include <armadillo>
#include <vector>

#include <xoos/types/float.h>

#include "segmentation/genomic-segments.h"
#include "sex.h"

namespace xoos::cnc {
using segmentation::GenomicSegment;
void FlagCallsByMeanMapq(std::vector<GenomicSegment>& segments, f64 mapq_cutoff_for_calls);
void FlagCallsByExpectedTotalCopyNumber(std::vector<GenomicSegment>& segments, const Sex& sex);
void FlagCallsByLength(std::vector<GenomicSegment>& segments, size_t cnv_length_flag_min_size);
void FlagChrYCallsIfFemale(std::vector<GenomicSegment>& segments);
}  // namespace xoos::cnc
