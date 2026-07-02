#pragma once
#include <optional>
#include <string>
#include <utility>

#include <taskflow/core/graph.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>

#include "baits.h"
#include "copy-number-caller/common/coverage-check.h"
#include "coverage.h"

namespace xoos::cnc {
class BamToCorrectedCoverageTaskFlowGraph {
 public:
  BamToCorrectedCoverageTaskFlowGraph(const fs::path& bam_file,
                                      const BaitRecords& baits,
                                      const std::optional<fs::path>& coverage_out_fname,
                                      const std::optional<fs::path>& corrected_coverage_out_fname,
                                      std::string exclude_flags,
                                      bool ignore_dn,
                                      f64 gc_correct_first_span,
                                      io::CommandLineInfo command_line_info)
      : _bam_file(bam_file),
        _baits(baits),
        _coverage_out_fname(coverage_out_fname),
        _corrected_coverage_out_fname(corrected_coverage_out_fname),
        _exclude_flags(std::move(exclude_flags)),
        _ignore_dn(ignore_dn),
        _gc_correct_first_span(gc_correct_first_span),
        _command_line_info(std::move(command_line_info)),
        _graph(BuildGraph()) {
  }

  // NOLINTNEXTLINE
  tf::Graph& graph() {
    return _graph;
  }

  CoverageRecords& GetResult() {
    return _res;
  }

  bool IsLowCoverage() const {
    return cnc::IsLowCoverage(_coverage_check_result);
  }

  const std::optional<CoverageCheckResult>& GetCoverageCheckResult() const {
    return _coverage_check_result;
  }

 private:
  tf::Graph BuildGraph();
  const fs::path& _bam_file;
  const BaitRecords& _baits;
  CoverageRecords _uncorrected_coverage;
  CoverageRecords _res;
  std::optional<fs::path> _coverage_out_fname;
  std::optional<fs::path> _corrected_coverage_out_fname;
  std::string _exclude_flags;
  bool _ignore_dn;
  f64 _gc_correct_first_span;
  io::CommandLineInfo _command_line_info;
  tf::Graph _graph;
  std::optional<CoverageCheckResult> _coverage_check_result;
};
}  // namespace xoos::cnc
