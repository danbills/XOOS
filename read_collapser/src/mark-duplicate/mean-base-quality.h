#pragma once

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/float.h>
#include <xoos/types/vec.h>

#include "io/alignment.h"

namespace xoos::read_collapser {

/**
 * Calculate the mean base quality of the bases in the alignment.
 */
f64 MeanBaseQ(const AlignmentPtr& alignment);

/**
 * Find the best representative alignment from @p alignments for duplicate marking.
 * The alignment with the highest mean base quality is selected. Full-length reads are preferred over partial reads.
 * Ties are broken by input order (first occurrence wins). Returns nullptr if @p alignments is empty.
 */
AlignmentPtr FindAlignmentWithMaxMeanBaseQ(const vec<AlignmentPtr>& alignments);

}  // namespace xoos::read_collapser
