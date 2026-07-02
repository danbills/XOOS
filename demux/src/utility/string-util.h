#pragma once

#include <string>
#include <string_view>

namespace xoos::demux {
std::string_view Substr(const std::string_view& str, std::string::size_type start_pos, std::string::size_type end_pos);

std::string Reverse(const std::string_view& str);
}  // namespace xoos::demux
