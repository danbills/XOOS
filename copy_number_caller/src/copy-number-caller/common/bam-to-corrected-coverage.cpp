#include "copy-number-caller/common/bam-to-corrected-coverage.h"

#include <taskflow/taskflow.hpp>

#include <xoos/gc_correct/gc-correct-taskflow-graph.h>
#include <xoos/log/logging.h>

#include "calculate-coverage/calculate-coverage-taskflow-graph.h"
#include "copy-number-caller/common/coverage-check.h"

namespace xoos::cnc {
tf::Graph BamToCorrectedCoverageTaskFlowGraph::BuildGraph() {
  tf::Graph graph;
  tf::FlowBuilder builder(graph);
  tf::Task coverage_task = builder.emplace([this](tf::Subflow& subflow) {
    CalculateCoverageTaskflowGraph calculate_coverage_tf_graph(_bam_file, _baits, _exclude_flags, _ignore_dn);
    subflow.composed_of(calculate_coverage_tf_graph);
    subflow.join();
    _uncorrected_coverage = calculate_coverage_tf_graph.GetResult();
    if (_coverage_out_fname.has_value()) {
      Logging::Info("Writing coverage to {}", _coverage_out_fname->string());
      std::ofstream ofs(_coverage_out_fname.value());
      _uncorrected_coverage.Write(ofs, _command_line_info);
    }
    // Check median coverage before proceeding to GC correction
    _coverage_check_result = CheckMedianCoverage(_uncorrected_coverage);
    if (!_coverage_check_result->is_sufficient) {
      Logging::Error(
          "Median coverage ({:.2f}x) is below minimum threshold ({:.0f}x). "
          "Skipping copy number calling.",
          _coverage_check_result->median_coverage,
          kMinMedianCoverageThreshold);
    }
  });
  tf::Task gc_correct_task = builder.emplace([this](tf::Subflow& subflow) {
    if (IsLowCoverage()) {
      Logging::Info("Skipping GC correction due to low coverage");
      return;
    }
    gc_correct::GCCorrectTaskFlowGraph gc_correct_taskflow_graph(_uncorrected_coverage.region,
                                                                 _uncorrected_coverage.count,
                                                                 _uncorrected_coverage.total_coverage,
                                                                 _baits.GetGCBias(),
                                                                 _baits.GetMappability(),
                                                                 _baits.GetOnTargetStatus(),
                                                                 _gc_correct_first_span);
    subflow.composed_of(gc_correct_taskflow_graph);
    subflow.join();
    gc_correct::GCCorrectResults& gc_correct_res = gc_correct_taskflow_graph.GetResult();
    _res.region = std::move(_uncorrected_coverage.region);
    _res.count = std::move(gc_correct_res.counts);
    _res.total_coverage = std::move(gc_correct_res.total_coverage);
    _res.average_coverage = std::move(gc_correct_res.average_coverage);
    _res.on_target = std::move(_uncorrected_coverage.on_target);
    _res.mean_mapping_quality = std::move(_uncorrected_coverage.mean_mapping_quality);
    _res.has_rescued_secondaries = _uncorrected_coverage.has_rescued_secondaries.load();
    _res.UpdateRegionMap();
    if (_corrected_coverage_out_fname.has_value()) {
      Logging::Info("Writing corrected coverage to {}", _corrected_coverage_out_fname->string());
      std::ofstream ofs(_corrected_coverage_out_fname.value());
      _res.Write(ofs, _command_line_info);
    }
  });
  coverage_task.precede(gc_correct_task);
  return graph;
}
}  // namespace xoos::cnc
