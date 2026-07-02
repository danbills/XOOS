#pragma once

#include <map>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

// This file contains template implementations of transparent functions that utilize the std::string key
// to enhance performance
// https://rules.sonarsource.com/cpp/tag/since-c++14/RSPEC-6045/

namespace xoos {

// Transparent hasher for string types. Enables heterogeneous lookup with
// std::string_view and const char* without constructing a temporary std::string.
struct StringHash {
  // NOLINTNEXTLINE(readability-identifier-naming) name is required by the C++ standard
  using is_transparent = void;

  size_t operator()(std::string_view sv) const noexcept;  // NOSONAR
};

// unordered containers with transparent functions
template <typename Tp, typename Alloc = std::allocator<Tp>>
using UnorderedSet = std::unordered_set<Tp, std::hash<Tp>, std::equal_to<>, Alloc>;
using StrUnorderedSet = std::unordered_set<std::string, StringHash, std::equal_to<>>;

template <typename Tp, typename Alloc = std::allocator<std::pair<const std::string, Tp>>>
using StrUnorderedMap = std::unordered_map<std::string, Tp, StringHash, std::equal_to<>, Alloc>;

// ordered containers with transparent functions
using StrSet = std::set<std::string, std::less<>>;

template <typename Tp, typename Alloc = std::allocator<std::pair<const std::string, Tp>>>
using StrMap = std::map<std::string, Tp, std::less<>, Alloc>;

}  // namespace xoos
