#pragma once

#include <IITree.h>

#include <string>
#include <unordered_map>
#include <vector>

#include "observations.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc::segmentation {
class IntervalTrees {
 public:
  IntervalTrees() = default;
  IntervalTrees(const IntervalTrees&) = default;
  IntervalTrees(IntervalTrees&&) = default;
  IntervalTrees& operator=(const IntervalTrees&) = default;
  IntervalTrees& operator=(IntervalTrees&&) = default;
  ~IntervalTrees() = default;
  explicit IntervalTrees(const std::vector<std::string>& regions);
  explicit IntervalTrees(const Observations& obvs);
  explicit IntervalTrees(const std::vector<GenomicSegment>& segments);
  std::vector<size_t> LookUp(std::string contig, size_t start, size_t end) const;

 private:
  // unordered map to store IITree (interval trees) for each contig
  // std::string: contig name
  // size_t index type for IITree
  // size_t data type for IITree
  using IntervalTreesT = std::unordered_map<std::string, IITree<size_t, size_t>>;
  IntervalTreesT _trees;
};
}  // namespace xoos::cnc::segmentation
