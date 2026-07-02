#include "segmentation/read-segments.h"

#include <algorithm>
#include <string>

#include <fmt/core.h>

#include <csv.hpp>

#include <xoos/enum/enum-util.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/util/string-functions.h>

#include "segmentation/genomic-segment-column-names.h"
#include "segmentation/segments-header.h"

namespace xoos::cnc::segmentation {

/**
 * @brief Extract the value of a CLI option from a serialized command line string.
 *
 * This parses the already-rendered command line stored in ##RocheCommandLine metadata.
 * CLI11 cannot be used here because we are operating on a flat string, not a live CLI11 app.
 * The format is always "--option value" with space separation (produced by RenderCli).
 */
static std::string ExtractCliOption(const std::string& command_line, const std::string& option_name) {
  auto tokens = xoos::string::Split(command_line, " ");
  auto it = std::find(tokens.begin(), tokens.end(), option_name);
  if (it == tokens.end() || ++it == tokens.end()) {
    return "";
  }
  return *it;
}

static std::vector<std::string> GetColumnNames(SegmentType mode) {
  switch (mode) {
    case SegmentType::kLogROnly:
      return kGenomicSegColsRequiredLogRSegments;
    case SegmentType::kBaf:
      return kGenomicSegColsRequiredBafSegments;
    case SegmentType::kGermlineLikelihood:
      return kGenomicSegColsRequiredGermline;
    case SegmentType::kSomaticWithBafLikelihood:
      return kGenomicSegColsRequiredSomaticAlleleSpecific;
    case SegmentType::kSomaticNoBafLikelihood:
      return kGenomicSegColsRequiredSomaticNoAlleleSpecific;
    case SegmentType::kSeed:
      return kGenomicSegColsRequiredSeed;
    default:
      Logging::Error("invalid SegmentType {}", enum_util::FormatEnumName(mode));
      throw std::runtime_error("invalid SegmentType");
  }
}

/** @brief
 *
 * we have to read a string here so that the parser can "magically" determine the format and know when to skip comments
 */
std::vector<GenomicSegment> ReadSegments(const fs::path& fname, SegmentType segment_type) {
  std::vector<GenomicSegment> ret;
  csv::CSVReader reader(fname.string());  //, format);
  auto col_names = GetColumnNames(segment_type);
  s32 row_num = 0;
  for (const csv::CSVRow& row : reader) {
    // check if all the required columns are present in this row
    auto row_col_names = row.get_col_names();
    GenomicSegment seg;
    for (const std::string& col : col_names) {
      if (std::find(row_col_names.begin(), row_col_names.end(), col) == row_col_names.end()) {
        Logging::Error("missing column {} in file {} at row number {}", col, fname.string(), row_num);
        throw std::runtime_error(fmt::format("missing column {}", col));
      }
      csv::CSVField val = row[col];
      seg.SetFieldOrThrowError(col, val.get<std::string>());
      if (segment_type == SegmentType::kGermlineLikelihood) {
        seg.purity = 0.99;
        seg.ploidy = 2;
      }
    }
    if (seg.end <= seg.start) {
      Logging::Error("SEG file error: end position must be greater than start position in file {} at row number {}",
                     fname.string(),
                     row_num);
      throw std::runtime_error("invalid segment");
    }
    ret.emplace_back(seg);
    row_num += 1;
  }
  Logging::Debug("read {} rows", row_num);
  return ret;
}

SegmentsHeader ReadHeaderFromSegments(const fs::path& fname) {
  SegmentsHeader ret;
  std::ifstream ifs(fname);
  std::string line;
  while (std::getline(ifs, line)) {
    if (line.find('#', 0) != 0) {
      break;
    }
    // New format: ##RocheCommandLine=<ID=...,Version="...",CommandLine="...">
    auto parsed = io::ParseCommandInfo(line);
    if (parsed.has_value()) {
      ret.program_name = parsed->name;
      ret.version = parsed->version;
      ret.command_line = parsed->command_line;
      if (ret.reference_file.empty()) {
        ret.reference_file = ExtractCliOption(parsed->command_line, "--reference");
      }
      continue;
    }
    // Old format: #KEY=VALUE
    if (line.find(kSegmentsHeaderProgramName, 1) == 1) {
      ret.program_name = line.substr(2 + kSegmentsHeaderProgramName.size());
    } else if (line.find(kSegmentsHeaderVersion, 1) == 1) {
      ret.version = line.substr(2 + kSegmentsHeaderVersion.size());
    } else if (line.find(kSegmentsHeaderReferenceFile, 1) == 1) {
      ret.reference_file = line.substr(2 + kSegmentsHeaderReferenceFile.size());
    } else if (line.find(kSegmentsHeaderCommandLine, 1) == 1) {
      ret.command_line = line.substr(2 + kSegmentsHeaderCommandLine.size());
    }
  }
  return ret;
}
}  // namespace xoos::cnc::segmentation
