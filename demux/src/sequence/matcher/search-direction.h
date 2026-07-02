#pragma once

namespace xoos::demux {

/**
 * Specifies the intended scanning direction for sequence matching.
 *
 * kForward  — match position relative to the left (5') end of the sequence.
 * kBackward — match position relative to the right (3') end of the sequence.
 *
 * The mechanism varies by implementation (e.g. pattern reversal, input iteration order).
 * See individual class documentation for details.
 */
enum class SearchDirection {
  kForward,
  kBackward,
};

}  // namespace xoos::demux
