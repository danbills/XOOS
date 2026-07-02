#pragma once
#include <string>

#include <taskflow/core/graph.hpp>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "baits.h"
#include "coverage.h"

namespace xoos::cnc {

class CalculateCoverageTaskflowGraph {
 public:
  CalculateCoverageTaskflowGraph(const fs::path& bam_file,
                                 const BaitRecords& baits,
                                 const std::string& exclude_flags,
                                 bool ignore_dn)
      : _bam_file(bam_file), _baits(baits), _exclude_flags(exclude_flags), _ignore_DN(ignore_dn), _graph(BuildGraph()) {
  }

  // NOLINTNEXTLINE
  tf::Graph& graph() {
    return _graph;
  }

  const CoverageRecords& GetResult() const {
    return _res;
  }

  CoverageRecords& GetResult() {
    return _res;
  }

 private:
  tf::Graph BuildGraph();
  const fs::path& _bam_file;
  const BaitRecords& _baits;
  const std::string& _exclude_flags;
  u32 _exclude_flags_int = 1796;
  bool _ignore_DN;
  s32 _step = 0;
  CoverageRecords _res;
  tf::Graph _graph;
};
}  // namespace xoos::cnc
