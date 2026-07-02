#include "calculate-coverage/calculate-coverage-taskflow-graph.h"

#include <htslib/sam.h>

#include <taskflow/algorithm/for_each.hpp>

#include <xoos/io/alignment-reader.h>
#include <xoos/log/logging.h>

#include "calculate-coverage/calculate-coverage.h"
#include "utility/utility-functions.h"

namespace xoos::cnc {
tf::Graph CalculateCoverageTaskflowGraph::BuildGraph() {
  tf::Graph graph;
  tf::FlowBuilder builder(graph);
  _res.region = _baits.GetRegions();
  _res.on_target = _baits.GetOnTargetStatus();
  size_t total_regions = _res.region.size();
  _res.total_coverage.zeros(total_regions);
  _res.count.zeros(total_regions);
  _res.mean_mapping_quality.zeros(total_regions);
  _res.has_duplicates = false;
  // _exclude_flags_int = (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);
  _exclude_flags_int = 0;
  s32 flags = bam_str2flag(_exclude_flags.c_str());
  if (flags < 0 || flags > ((BAM_FSUPPLEMENTARY << 1) - 1)) {
    throw std::runtime_error("Error: Flag value not supported");
  }
  _exclude_flags_int |= flags;
  auto region_tuple = ParseRegionString(_res.region[0]);
  auto start = static_cast<s32>(std::get<1>(region_tuple));
  auto end = static_cast<s32>(std::get<2>(region_tuple));
  s32 region_size = end - start;
  // step size is optimized for WGS, which is around 1 - 2 million bp per thread
  _step = 1500000 / region_size;
  if (_step > static_cast<s32>(total_regions)) {
    _step = static_cast<s32>(total_regions);
  }
  if (_step < 1) {
    _step = 1;
  }
  Logging::Debug("Step size {}, region size {}", _step, region_size);
  Logging::Debug("Exclude flags {} {}", _exclude_flags, _exclude_flags_int);
  tf::Task tf = builder.for_each_index(0, static_cast<s32>(total_regions), _step, [this](s32 i) {
    const io::AlignmentReader reader = io::OpenAlignmentReader(_bam_file);
    CalculateCoverageRegion(reader, _res, _exclude_flags_int, _ignore_DN, i, _step);
  });
  tf::Task empty_task = builder.emplace([]() {});
  tf.precede(empty_task);
  return graph;
}
}  // namespace xoos::cnc
