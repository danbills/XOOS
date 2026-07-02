#include "utility/utility-functions.h"

#include <fstream>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/util/parse-int.h>

namespace xoos::cnc {

bool IsTrueString(const std::string& s) {
  if (s == "1" || s == "TRUE" || s == "true") {
    return true;
  } else if (s == "0" || s == "FALSE" || s == "false") {
    return false;
  } else {
    throw std::runtime_error("Invalid boolean string: " + s);
  }
}

size_t StringToNonNegativeUIntOrThrow(const std::string& s) {
  return util::ParseU64(s);
}

/**
 * @brief parses region string into contig, start and end. start is 0-based and end is 1-based exclusive (as in a BED
 * file)
 * @param r  region string
 * @return tuple - contig, start and end
 */
std::tuple<std::string, size_t, size_t> ParseRegionString(const std::string& r) {
  auto colon_pos = r.find(':');
  auto dash_pos = r.find('-');
  auto contig = r.substr(0, colon_pos);
  auto start = util::ParseU64(r.substr(colon_pos + 1, dash_pos - colon_pos - 1));
  if (start == 0) {
    throw error::Error("start position in region string cannot be 0");
  }
  start -= 1;
  auto end = util::ParseU64(r.substr(dash_pos + 1));
  return std::make_tuple(contig, start, end);
}

bool IsEqualContig(std::string_view region, std::string_view contig) {
  const auto colon_pos = region.find(':');
  const auto prefix_len = (colon_pos == std::string_view::npos) ? region.size() : colon_pos;
  return prefix_len == contig.size() && region.compare(0, prefix_len, contig) == 0;
}

std::vector<fs::path> PathsFromFile(const fs::path& fname) {
  std::ifstream ifs(fname);
  std::vector<fs::path> paths;
  std::string line;
  while (std::getline(ifs, line)) {
    paths.emplace_back(line);
  }
  return paths;
}

}  // namespace xoos::cnc
