#include "segmentation/interval-trees.h"

#include <algorithm>

#include <xoos/log/logging.h>

#include "utility/utility-functions.h"

namespace xoos::cnc::segmentation {

/**
 * @brief creates interval trees for list of region strigns
 * @param regions region strings of the format contig:start-end (1-based closed)
 * @return map, key=contig, item=interval tree for segments in the contig. data value for each node in the tree
 * corresponds to index w/i `segments`
 */
IntervalTrees::IntervalTrees(const std::vector<std::string>& regions) {
  for (const auto& region : regions) {
    const auto& [contig, start, end] = ParseRegionString(region);
    _trees[contig].add(start, end, regions.size());  // ParseRegionString already converts to 0-based half-open
  }
  for (auto& [contig_name, interval_tree] : _trees) {
    interval_tree.index();
  }
}

/**
 * @brief creates interval trees for Observations
 * @param obvs input observations
 * @return map, key=contig, item=interval tree for segments in the contig. data value for each node in the tree
 * corresponds to index w/i `segments`
 */
IntervalTrees::IntervalTrees(const Observations& obvs) {
  for (size_t i = 0; i < obvs.contigs.size(); ++i) {
    _trees[obvs.contigs[i]].add(obvs.starts[i], obvs.ends[i], i);
  }
  for (auto& [contig_name, interval_tree] : _trees) {
    interval_tree.index();
  }
}

/**
 * @brief creates interval trees for genomic segments
 * @param segments input segments
 * @return map, key=contig, item=interval tree for segments in the contigof corresponding nod. data value for each node
 * in the tree corresponds to index w/i `segments`
 */
IntervalTrees::IntervalTrees(const std::vector<GenomicSegment>& segments) {
  for (size_t i = 0; i < segments.size(); ++i) {
    _trees[segments[i].contig].add(segments[i].start, segments[i].end, i);
  }
  for (auto& pair : _trees) {
    pair.second.index();
  }
}

/**
 * @brief  looks up a <contig start end> pair inside _trees.
 * @param contig
 * @param start
 * @param end
 * @return idxs within original from which tree was constructed that represent overlaps to input
 */
std::vector<size_t> IntervalTrees::LookUp(std::string contig, size_t start, size_t end) const {
  const auto tree_ptr = _trees.find(contig);
  if (tree_ptr == _trees.end()) {
    Logging::Debug("contig {} is not observed in IntervalTrees", contig);
    return {};
  }
  const auto& tree = tree_ptr->second;
  // overlaps will contain the indexes in seed_segments that obvs[i] overlaps
  std::vector<size_t> overlaps;
  // IITree overlap() method fills in `overlaps` std::vector instead of returning a value
  tree.overlap(start, end, overlaps);
  // TODO: make sure segments span whole chromomosome so that this never happens?
  if (overlaps.empty()) {
    Logging::Debug("observation {}:{}-{} does not appear in IntervalTrees. skipping...", contig, start + 1, end);
    return {};
  }
  std::vector<size_t> idxs(overlaps.size());
  for (size_t i = 0; i < overlaps.size(); ++i) {
    idxs[i] = tree.data(overlaps[i]);
  }
  std::sort(idxs.begin(), idxs.end());
  return idxs;
}
}  // namespace xoos::cnc::segmentation
