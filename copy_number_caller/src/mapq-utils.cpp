#include "mapq-utils.h"

#include <xoos/error/error.h>

#include "segmentation/interval-trees.h"
#include "utility/utility-functions.h"

namespace xoos::cnc {
using segmentation::IntervalTrees;

/**
 * @brief given coverage records with mean mapping qualities, return observations with those mapping qualities
 * @param coverage object with mean_mapping_quality field populated
 * @return Observations object where the observations consist of the mean_mapping_quality field from `coverage`
 */
Observations GetAvgMapqsFromCoverage(const CoverageRecords& coverage) {
  Observations mapqs;
  mapqs.contigs.resize(coverage.region.size());
  mapqs.starts.resize(coverage.region.size());
  mapqs.ends.resize(coverage.region.size());
  mapqs.obvs.resize(coverage.region.size());
  for (size_t i = 0; i < coverage.region.size(); ++i) {
    auto [contig, start, end] = ParseRegionString(coverage.region[i]);
    mapqs.contigs[i] = contig;
    mapqs.starts[i] = start;
    mapqs.ends[i] = end;
    mapqs.obvs[i] = coverage.mean_mapping_quality[i];
  }
  return mapqs;
}

/**
 * @brief given a set of segments (each segment represents a collection of observations), and a list of MAPQs per
 * observation, return a list of mean MAPQs per segment
 * @param mapqs
 * @param segments
 * @return
 */
arma::vec GetAvgMeanMapqPerSegment(const std::vector<GenomicSegment>& segments, const Observations& mapqs) {
  arma::vec avg_mean_mapqs(segments.size());
  IntervalTrees mapqs_trees(mapqs);
  size_t i = 0;
  for (const GenomicSegment& seg : segments) {
    if (!seg.num_obs.has_value() || seg.num_obs.value() == 0) {
      throw error::Error(
          "Segment {}:{}-{} must have non-zero num_obs to calculate mean MAPQ", seg.contig, seg.start, seg.end);
    }
    std::vector<size_t> idxs = mapqs_trees.LookUp(seg.contig, seg.start, seg.end);
    f64 avg_mean_mapq = 0;
    for (size_t idx : idxs) {
      avg_mean_mapq += mapqs.obvs[idx];
    }
    avg_mean_mapq = avg_mean_mapq / static_cast<f64>(seg.num_obs.value());
    avg_mean_mapqs[i++] = avg_mean_mapq;
  }
  return avg_mean_mapqs;
}

/**
 * @brief given a set of segments (each segment represents a collection of observations), and a list of MAPQs per
 * observation, assign each segment's avg_mean_mapq field with the average mean mapq of all the targets in that segment
 * @param mapqs, non-const
 * @param segments, const
 * @return
 */
std::vector<GenomicSegment>& AssignAvgMeanMapqPerSegment(std::vector<GenomicSegment>& segments,
                                                         const Observations& mapqs) {
  IntervalTrees mapqs_trees(mapqs);
  for (GenomicSegment& seg : segments) {
    if (!seg.num_obs.has_value() || seg.num_obs.value() == 0) {
      throw error::Error(
          "Segment {}:{}-{} must have non-zero num_obs to calculate mean MAPQ", seg.contig, seg.start, seg.end);
    }
    std::vector<size_t> idxs = mapqs_trees.LookUp(seg.contig, seg.start, seg.end);
    f64 avg_mean_mapq = 0;
    f64 n = 0;
    for (size_t idx : idxs) {
      if (!std::isnan(mapqs.obvs[idx])) {
        avg_mean_mapq += mapqs.obvs[idx];
        n += 1;
      }
    }
    if (n > 0) {
      seg.avg_mean_mapq = avg_mean_mapq / n;
    } else {
      seg.avg_mean_mapq = std::nullopt;
    }
  }
  return segments;
}

}  // namespace xoos::cnc
