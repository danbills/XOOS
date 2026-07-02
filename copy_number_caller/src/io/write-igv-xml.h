#pragma once
#include <optional>

#include <xoos/types/fs.h>

namespace xoos::cnc {
void WriteIGVXML(const fs::path& output_xml_path,
                 const std::optional<fs::path>& logr_file,
                 const std::optional<fs::path>& bigwig_file);
}  // namespace xoos::cnc
