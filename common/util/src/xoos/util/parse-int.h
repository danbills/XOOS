#pragma once
#include <concepts>
#include <string_view>

#include <xoos/types/int.h>

namespace xoos::util {

/// Parse a string as an unsigned integer of type T.
/// Throws on empty input, non-digit characters, or overflow beyond T's range.
template <std::unsigned_integral T>
T ParseUnsigned(std::string_view value);

/// Parse a string as a signed integer of type T.
/// Accepts an optional leading minus sign.
/// Throws on empty input, invalid characters, or overflow beyond T's range.
template <std::signed_integral T>
T ParseSigned(std::string_view value);

// Convenience aliases preserving the original API.

inline u64 ParseU64(const std::string_view value) {
  return ParseUnsigned<u64>(value);
}

inline u32 ParseU32(const std::string_view value) {
  return ParseUnsigned<u32>(value);
}

inline u16 ParseU16(const std::string_view value) {
  return ParseUnsigned<u16>(value);
}

inline u8 ParseU8(const std::string_view value) {
  return ParseUnsigned<u8>(value);
}

inline s64 ParseS64(const std::string_view value) {
  return ParseSigned<s64>(value);
}

inline s32 ParseS32(const std::string_view value) {
  return ParseSigned<s32>(value);
}

inline s16 ParseS16(const std::string_view value) {
  return ParseSigned<s16>(value);
}

inline s8 ParseS8(const std::string_view value) {
  return ParseSigned<s8>(value);
}

}  // namespace xoos::util
