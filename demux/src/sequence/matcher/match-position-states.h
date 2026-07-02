#pragma once

#include <xoos/types/int.h>

#include <limits>

namespace xoos::demux {

constexpr s32 kNoMatchPosition = -1;
constexpr s32 kPositionAfterMatchAdjustment = 1;
constexpr u64 kEndPosition = std::numeric_limits<u64>::max();

}  // namespace xoos::demux
