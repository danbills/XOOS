#include "xoos/util/string-functions.h"

#include <cstdint>
#include <vector>

namespace xoos::string {

std::vector<std::string> Split(const std::string& s, const std::string& delimiter) {
  size_t pos_start = 0;
  size_t pos_end;
  const size_t delim_len = delimiter.length();
  std::string token;
  std::vector<std::string> res;
  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
    token = s.substr(pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res.push_back(token);
  }
  res.push_back(s.substr(pos_start));
  return res;
}

std::string Join(const std::vector<std::string>& input_vector, const std::string& delimiter) {
  std::string result;
  uint32_t counter = 0;
  for (const auto& element : input_vector) {
    if (counter == input_vector.size() - 1) {
      result += element;
    } else {
      result += element + delimiter;
    }
    ++counter;
  }
  return result;
}

}  // namespace xoos::string
