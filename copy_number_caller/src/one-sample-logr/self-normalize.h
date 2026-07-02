#pragma once

#include "coverage.h"
#include "observations.h"

namespace xoos::cnc {
Observations SelfNormalizeCounts(const CoverageRecords& coverage);
}
