#pragma once
#include <string_view>

namespace xoos::demux {

constexpr std::string_view kNameKey = "name";
constexpr std::string_view kTypeKey = "type";
constexpr std::string_view kAdapter5pKey = "adapter_5p";
constexpr std::string_view kAdapter3pKey = "adapter_3p";
constexpr std::string_view kBait5pKey = "bait_5p";
constexpr std::string_view kBait3pKey = "bait_3p";

constexpr std::string_view kSequenceKey = "sequence";
constexpr std::string_view kSequencesKey = "sequences";
constexpr std::string_view kMaxEditDistanceKey = "max_edit_distance";
constexpr std::string_view kTransformKey = "transform";
constexpr std::string_view kPrefixKey = "prefix";
constexpr std::string_view kSuffixKey = "suffix";

}  // namespace xoos::demux
