#pragma once

#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "io/alignment.h"

namespace xoos::read_collapser {

struct ClusterId;

void DownsampleReadsInCluster(vec<AlignmentPtr>& reads, u32 max_reads, const ClusterId& cluster_id);

}  // namespace xoos::read_collapser
