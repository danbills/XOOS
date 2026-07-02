#pragma once
#include <armadillo>
#include <vector>

#include "coverage.h"
#include "observations.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/interval-trees.h"

namespace xoos::cnc {
using segmentation::GenomicSegment;
using segmentation::IntervalTrees;

static const std::string kMeanMapqColumnName = "MeanMapQ";

Observations GetAvgMapqsFromCoverage(const CoverageRecords& coverage);
arma::vec GetAvgMeanMapqPerSegment(const std::vector<GenomicSegment>& segments, const Observations& mapqs);
std::vector<GenomicSegment>& AssignAvgMeanMapqPerSegment(std::vector<GenomicSegment>& segments,
                                                         const Observations& mapqs);

}  // namespace xoos::cnc
