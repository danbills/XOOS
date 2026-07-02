#include "likelihood/likelihood-flags.h"

#include <ranges>

#include <xoos/enum/enum-util.h>

namespace xoos::cnc {
std::vector<std::string> StringToFlags(const std::string& str) {
  std::vector<std::string> ret;
  for (const auto& word : std::views::split(str, std::string_view{";"})) {
    std::string flag;
    std::ranges::copy(word, std::back_inserter(flag));
    if (std::find(kLikelihoodFlagsAll.begin(), kLikelihoodFlagsAll.end(), flag) == kLikelihoodFlagsAll.end()) {
      Logging::Error("invalid flag: {}!", flag);
      throw std::runtime_error("flag parsing error");
    } else if (flag != kLikelihoodFlagNone) {
      ret.emplace_back(flag);
    }
  }
  return ret;
}

std::string FlagsToString(const std::vector<std::string>& flags) {
  if (flags.empty()) {
    return kLikelihoodFlagNone;
  }
  std::string flags_str{};
  for (const auto& flag : flags) {
    if (flag != kLikelihoodFlagNone) {
      flags_str += flag + ";";
    }
  }
  flags_str = flags_str.substr(0, flags_str.size() - 1);  // remove the last comma
  return flags_str;
}

}  // namespace xoos::cnc
