#include "consensus/base-encoder.h"

#include <xoos/log/logging.h>

namespace xoos::read_collapser {

BaseIndex ToIndex(const char base) {
  using enum BaseIndex;
  switch (base) {
    case kBaseA:
      return kIndexA;
    case kBaseC:
      return kIndexC;
    case kBaseG:
      return kIndexG;
    case kBaseT:
      return kIndexT;
    case kBaseGap:
      return kIndexGap;
    // N and P are valid consensus matrix symbols but callers (majority-voting.cpp)
    // filter them out before calling ToIndex. Mapped here for contract completeness.
    case kBaseN:
      return kIndexN;
    case kBaseP:
      return kIndexP;
    default: {
      XOOS_LOG_WARN_ONCE("ToIndex: unexpected base character '{}'", base);
      return kIndexGap;
    }
  }
}

size_t ToSizeTIndex(const char base) {
  return static_cast<size_t>(ToIndex(base));
}

char ToBase(const BaseIndex index) {
  using enum BaseIndex;
  switch (index) {
    case kIndexA:
      return kBaseA;
    case kIndexC:
      return kBaseC;
    case kIndexG:
      return kBaseG;
    case kIndexT:
      return kBaseT;
    case kIndexGap:
      return kBaseGap;
    // N and P are valid consensus matrix symbols but callers (majority-voting.cpp)
    // filter them out before calling ToBase. Mapped here for contract completeness.
    case kIndexN:
      return kBaseN;
    case kIndexP:
      return kBaseP;
    default: {
      XOOS_LOG_WARN_ONCE("ToBase: unexpected BaseIndex value {}", static_cast<size_t>(index));
      return kBaseGap;
    }
  }
}

}  // namespace xoos::read_collapser
