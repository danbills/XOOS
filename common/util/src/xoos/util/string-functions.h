#pragma once

#include <algorithm>
#include <locale>
#include <string>
#include <vector>

namespace xoos::string {

// string splitting function missing in standard library
std::vector<std::string> Split(const std::string& s, const std::string& delimiter);

// trim from start (in place)
static inline void LeftTrim(std::string& s) {
  s.erase(s.begin(),
          std::find_if(s.begin(), s.end(), [](const char ch) { return !std::isspace(ch, std::locale::classic()); }));
}

// trim from end (in place)
static inline void RightTrim(std::string& s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](const char ch) { return !std::isspace(ch, std::locale::classic()); })
              .base(),
          s.end());
}

// trim from both ends (in place)
static inline void Trim(std::string& s) {
  RightTrim(s);
  LeftTrim(s);
}

// trim from start (copying)
static inline std::string LeftTrimCopy(std::string s) {
  LeftTrim(s);
  return s;
}

// trim from end (copying)
static inline std::string RightTrimCopy(std::string s) {
  RightTrim(s);
  return s;
}

// trim from both ends (copying)
static inline std::string TrimCopy(std::string s) {
  Trim(s);
  return s;
}

std::string Join(const std::vector<std::string>& input_vector, const std::string& delimiter);

static inline bool EndsWithOneOf(const std::string& name, const std::vector<std::string>& suffixes) {
  std::string ps{name};
  std::ranges::transform(ps, std::begin(ps), [](const char c) { return std::tolower(c, std::locale::classic()); });
  return std::ranges::any_of(suffixes, [&ps](const auto& suffix) { return ps.ends_with(suffix); });
}

}  // namespace xoos::string
