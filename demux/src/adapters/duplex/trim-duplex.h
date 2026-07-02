#pragma once

#include <xoos/types/int.h>

#include <string_view>

#include "io/read-record.h"
#include "sequence/matcher/bitap.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

struct DuplexMetrics;

// TODO: hardcoded values should be removed, either taken from configuration or dynamically determined
// this should be based on the length of the sequence being searched for
constexpr s32 kBitapErrors = 4;

class TrimDuplex {
 public:
  explicit TrimDuplex(SeqMatcher runway_5p_matcher, SeqMatcher start_matcher);

  ~TrimDuplex() = default;

  // Find the start adapter in the consensus read
  s32 FindStartAdapterInConsensus(const FixedReadRecord& record) const;

  std::string_view Runway5pSequence() const;

 private:
  SeqMatcher _runway_5p_matcher;
  SeqMatcher _start_matcher;
  // Encapsulation of Bitap alignment algorithm to find start and end markers.
  // Edit distance tolerances are set to at least an error rate > 20% and at most 30% (~25% is target).
  // TODO: consider making this based on probability of random matches rather than fixed error rates (requires testing)
  // 4/18 = 22% error rate as 18bp is the expect length of the start adapter
  const Bitap<4> _kStart;
  // 3/12 = 25% error rate but 12 is a length determined empirically to support truncated adapters
  // Only searched for if the full-length start adapter is not found.
  const Bitap<3> _kShortStart;

  // Helper struct to pass around some results of the duplex trim.
  struct TrimResults {
    u32 length{0};
  };
};
}  // namespace xoos::demux
