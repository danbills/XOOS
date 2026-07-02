#include "sex.h"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <vector>

#include <fmt/core.h>

#include <xoos/log/logging.h>
#include <xoos/types/vec.h>

#include "utility/utility-functions.h"

namespace xoos::cnc {
// Support both hg38 (chr-prefixed) and hg19 (unprefixed) contig names.
const vec<std::string> kAllosomes({"chrX", "chrY", "X", "Y"});
const vec<std::string> kChromYNames({"chrY", "Y"});

bool SkipBasedOnSexChromHandling(const std::string& region, SexChromHandling sex_handling) {
  switch (sex_handling) {
    case SexChromHandling::kNoSexChrom: {
      return IsInAllosome(region);
    }
    case SexChromHandling::kOnlySexChrom: {
      return !IsInAllosome(region);
    }
    case SexChromHandling::kNoChromY: {
      return IsInChromY(region);
    }
    default: {
      return false;
    }
  }
  return false;
}

/**
 * @brief checks if a region string describes a region within a sex chromosome.
 * Uses exact contig matching to avoid false positives from alt contigs
 * (e.g., chrX_random, chrXY). Supports both hg38 (chrX/chrY) and hg19 (X/Y).
 * @param r  region string (e.g., "chrX:1-100" or "X:1-100")
 * @return true if the contig is exactly chrX, chrY, X, or Y
 */
bool IsInAllosome(std::string_view r) {
  return std::ranges::any_of(kAllosomes, [r](const std::string_view allosome) { return IsEqualContig(r, allosome); });
}

/**
 * @brief checks if a region string describes a region in chromosome Y.
 * Uses exact contig matching. Supports both hg38 (chrY) and hg19 (Y).
 * @param r  region string
 * @return true if the contig is exactly chrY or Y
 */
bool IsInChromY(std::string_view r) {
  return std::ranges::any_of(kChromYNames, [r](const std::string_view name) { return IsEqualContig(r, name); });
}

const std::vector<char> kPossibleSexChars{'M', 'F', 'N'};
const std::vector<std::string> kPossibleSexStrings{"M", "F", "N"};
const std::map<std::string, Sex> kCharToSex{{"M", Sex::kMale}, {"F", Sex::kFemale}, {"N", Sex::kUnknown}};
const std::map<Sex, std::string> kSexToChar{{Sex::kMale, "M"}, {Sex::kFemale, "F"}, {Sex::kUnknown, "N"}};

std::string SexToStr(Sex sex) {
  std::map<Sex, std::string> sex_to_str{{Sex::kMale, "M"}, {Sex::kFemale, "F"}, {Sex::kUnknown, "N"}};
  return sex_to_str[sex];
}

Sex ParseSexOption(std::string sex) {
  if (sex.empty()) {
    Logging::Info("no sex file or character specified, defaulting to 'N' (no likelihood model changes based on sex)");
    return Sex::kUnknown;
  }
  if (sex.size() != 1) {
    throw std::runtime_error("invalid file containing sex character, or invalid value for sex");
  }
  auto it = kCharToSex.find(sex.substr(0, 1));
  if (it == kCharToSex.end()) {
    throw std::runtime_error(fmt::format("invalid sex character: {}", sex[0]));
  }
  return it->second;
}
}  // namespace xoos::cnc
