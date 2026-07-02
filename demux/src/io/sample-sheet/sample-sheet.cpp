#include "sample-sheet.h"

#include <xoos/error/error.h>
#include <xoos/types/str-container.h>

#include <csv.hpp>

#include "task/flow-context.h"

namespace xoos::demux {

SidSeqNamePairs ReadSampleSheet(const fs::path& path) {
  SidSeqNamePairs sid_seq_name_pairs;
  StrUnorderedMap<std::string> sid2sample_name;
  StrUnorderedMap<std::string> sample_name2sid;
  csv::CSVFormat format;
  format.delimiter({',', '\t'});
  auto reader = csv::CSVReader(path.string(), format);
  for (const auto& row : reader) {
    auto sample_name = row[kSampleNameColumnHeader].get<std::string>();
    auto sample_sid = row[kSampleSidColumnHeader].get<std::string>();

    if (sample_name == kReservedRawFailedName) {
      throw error::Error("Sample name '{}' is reserved for unassigned reads and cannot be used in the sample sheet",
                         kReservedRawFailedName);
    }

    const auto it = sid2sample_name.find(sample_sid);
    if (it != sid2sample_name.end()) {
      throw error::Error("Found duplicate SID in sample sheet for '{}' and '{}'", sample_name, it->second);
    }
    sid2sample_name[sample_sid] = sample_name;

    const auto it_name = sample_name2sid.find(sample_name);
    if (it_name != sample_name2sid.end()) {
      throw error::Error("Found duplicate sample name in sample sheet for '{}' and '{}'", sample_sid, it_name->second);
    }
    sample_name2sid[sample_name] = sample_sid;

    sid_seq_name_pairs.emplace_back(sample_sid, sample_name);
  }
  return sid_seq_name_pairs;
}

BarcodePool LoadSampleSheet(const fs::path& sample_sheet) {
  auto result = BarcodePool{};
  const auto entries = ReadSampleSheet(sample_sheet);

  u32 id = 0;
  for (const auto& [sequence, name] : entries) {
    result.emplace_back(id, sequence, name);
    ++id;
  }
  return result;
}

}  // namespace xoos::demux
