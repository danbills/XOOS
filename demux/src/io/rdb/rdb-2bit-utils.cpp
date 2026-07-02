#include "rdb-2bit-utils.h"

#include <xoos/util/sequence-functions.h>

#include <algorithm>
#include <string>
#include <vector>

#include "sequence/encoding-constants.h"

namespace xoos::demux {

constexpr std::vector<u8> GenerateReverseLookupTable() {
  std::vector<u8> lookup_table(256);

  for (int i = 0; i < 256; ++i) {
    lookup_table[i] = (i & 0b11) << 6 | (i & 0b1100) << 2 | (i & 0b110000) >> 2 | (i & 0b11000000) >> 6;
  }

  return lookup_table;
}

static const std::vector<u8> kReverseLookupTable = GenerateReverseLookupTable();

void Reverse2BitOrder(u8* sequence, std::size_t sequence_length) {
  for (std::size_t i = 0; i < sequence_length; ++i) {
    sequence[i] = kReverseLookupTable[sequence[i]];
  }
}

constexpr std::vector<std::string> GenerateDnaDecodeTable() {
  std::vector<std::string> lookup_table(256);
  for (std::size_t i = 0; i < lookup_table.size(); ++i) {
    std::string decoded_bases;
    for (u32 j = 0; j < kBasesPerByte; ++j) {
      decoded_bases += sequence::kDnaAlphabet[(i >> (kBitsPerEncodedBase * j)) & kBaseMask];
    }
    lookup_table[i] = decoded_bases;
  }
  return lookup_table;
}

static const std::vector<std::string> kDnaLookupTable = GenerateDnaDecodeTable();

void DecodeDnaBases(const u8* const encoded, char* const buffer, const std::size_t length) {
  char* out = buffer;
  for (std::size_t i = 0; i < length >> kBitsPerEncodedBase; ++i) {
    const auto byte = encoded[i];
    const std::string& decoded_bases = kDnaLookupTable[byte];
    out = std::copy_n(decoded_bases.c_str(), decoded_bases.size(), out);
  }
  const std::size_t remaining = length & kBaseMask;
  if (remaining != 0) {
    const auto byte = encoded[length >> kBitsPerEncodedBase];
    const std::string& decoded_bases = kDnaLookupTable[byte];
    std::copy_n(decoded_bases.c_str(), remaining, out);
  }
}

constexpr std::vector<std::string> GenerateQualDecodeTable() {
  constexpr char kQual[] = {'+', '5', '?', 'I'};

  std::vector<std::string> lookup_table(256);
  for (int i = 0; i < 256; ++i) {
    std::string decoded_quals;
    for (s32 j = kBasesPerByte - 1; j >= 0; --j) {
      decoded_quals += kQual[(i >> (kBitsPerEncodedBase * static_cast<u32>(j))) & kBaseMask];
    }
    lookup_table[i] = decoded_quals;
  }
  return lookup_table;
}

static const std::vector<std::string> kQualLookupTable = GenerateQualDecodeTable();

void DecodeQualScores(const u8* const encoded, char* const buffer, const std::size_t length) {
  char* out = buffer;
  for (std::size_t i = 0; i < length >> kBitsPerEncodedBase; ++i) {
    uint8_t byte = encoded[i];
    const std::string& decoded_qual = kQualLookupTable[byte];
    out = std::copy_n(decoded_qual.c_str(), decoded_qual.size(), out);
  }
  const auto remaining = length & kBaseMask;
  if (remaining != 0) {
    const auto byte = encoded[length >> kBitsPerEncodedBase];
    const std::string& decoded_qual = kQualLookupTable[byte];
    std::copy_n(decoded_qual.c_str(), remaining, out);
  }
}

}  // namespace xoos::demux
