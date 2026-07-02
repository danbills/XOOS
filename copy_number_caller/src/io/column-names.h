#pragma once
#include <string>

#include <xoos/io/metadata-util.h>
#include <xoos/types/vec.h>

namespace xoos::cnc {
const std::string kColumnContig = "Contig";
const std::string kColumnContigAsComment = "#Contig";
const std::string kColumnStart = "Start";
const std::string kColumnEnd = "End";
const std::string kColumnMeanMapq = "MeanMapQ";
const std::string kColumnTotalCoverage = "TotalCoverage";
const std::string kColumnCounts = "Counts";
const std::string kColumnOnTarget = "OnTarget";
const std::string kColumnGCBias = "GCBias";
const std::string kColumnMappability = "Mappability";
const std::string kColumnLogRatio = "LogR";
const std::string kColumnRefAD = "Ref_AD";
const std::string kColumnAltAD = "Alt_AD";
const std::string kColumnBAF = "BAF";
}  // namespace xoos::cnc
