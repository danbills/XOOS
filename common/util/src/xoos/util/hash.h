#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <string_view>

namespace xoos::util::hash {

// Below hash functions are very similar to:
// D. Eastlake, “The FNV Non-Cryptographic Hash Algorithm,” Internet-Draft draft-eastlake-fnv-31,
// Internet Engineering Task Force, 2024. [Online]. Available:
// https://www.ietf.org/archive/id/draft-eastlake-fnv-31.html.

template <typename... Ts>
std::uint64_t HashCombine(Ts... values) {
  constexpr std::uint64_t kPrime{0x100000001b3};
  std::uint64_t result{0xcbf29ce484222325};
  for (const auto& value : {values...}) {
    result = (result * kPrime) ^ value;
  }
  return result;
}

template <typename... Ts>
std::uint64_t Hash(Ts... values) {
  constexpr std::uint64_t kPrime{0x100000001b3};
  std::uint64_t result{0xcbf29ce484222325};
  ([&] { result = (result * kPrime) ^ std::hash<Ts>{}(values); }(), ...);  // NOLINT
  return result;
}

template <typename Iter>
std::uint64_t HashRange(Iter begin, Iter end) {
  using ValueType = typename std::iterator_traits<Iter>::value_type;
  constexpr std::uint64_t kPrime{0x100000001b3};
  std::uint64_t result{0xcbf29ce484222325};
  for (auto iter = begin; iter != end; ++iter) {
    result = (result * kPrime) ^ std::hash<ValueType>{}(*iter);
  }
  return result;
}

/**
 * Transparent hash function for strings that supports heterogeneous lookup in unordered containers.
 * This allows us to use std::string_view or const char* as keys for lookup in an unordered_set of
 * std::string without needing to construct temporary std::string objects.
 */
struct TransparentStringHash {
  // NOLINTNEXTLINE(readability-identifier-naming)
  using is_transparent [[maybe_unused]] = void;

  [[nodiscard]] std::size_t operator()(const std::string& str) const {
    return std::hash<std::string>{}(str);
  }

  [[nodiscard]] std::size_t operator()(std::string_view str) const {
    return std::hash<std::string_view>{}(str);
  }

  [[nodiscard]] std::size_t operator()(const char* str) const {
    return std::hash<std::string_view>{}(str);
  }
};

}  // namespace xoos::util::hash
