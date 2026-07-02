#include "xoos/types/str-container.h"

#include <functional>

namespace xoos {

size_t StringHash::operator()(const std::string_view sv) const noexcept {
  return std::hash<std::string_view>{}(sv);
}

}  // namespace xoos
