#pragma once
#include <string>

namespace xoos::cnc::segmentation {
const std::string kSegmentsHeaderProgramName = "PROGRAM_NAME";
const std::string kSegmentsHeaderVersion = "VERSION";
const std::string kSegmentsHeaderReferenceFile = "REFERENCE_FILE";
const std::string kSegmentsHeaderCommandLine = "COMMAND_LINE";

struct SegmentsHeader {
  std::string program_name;
  std::string version;
  std::string reference_file;
  std::string command_line;
};
}  // namespace xoos::cnc::segmentation
