#pragma once

#include <map>
#include <string>
#include <string_view>

namespace xoos::cnc {

enum class SexChromHandling {
  kAllChrom,
  kNoSexChrom,
  kOnlySexChrom,
  kNoChromY
};
enum class Sex {
  kMale,
  kFemale,
  kUnknown
};
extern const std::map<std::string, Sex> kCharToSex;
extern const std::map<Sex, std::string> kSexToChar;
bool IsInAllosome(std::string_view r);
bool IsInChromY(std::string_view r);
bool SkipBasedOnSexChromHandling(const std::string& region, SexChromHandling sex_handling);
std::string SexToStr(Sex sex);
Sex ParseSexOption(std::string sex);

}  // namespace xoos::cnc
