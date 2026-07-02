#include "task/demux.h"

#include <xoos/log/logging.h>

#include "adapters/duplex/trim-info-duplex.h"
#include "metrics/duplex-metrics.h"
#include "metrics/read-completeness.h"
#include "metrics/simplex-metrics.h"
#include "task/flow-context.h"
#include "utility/stop-watch.h"

namespace xoos::demux {

Demux::Demux(FlowContext& exec, size_t batch_nr) : Task(exec, batch_nr, fmt::format("Demux{}", batch_nr)) {}

/**
 * @brief RunSimplexAdapter is a helper function to run the demuxing and trimming for simplex adapters
 * @param read the read to be demultiplexed and trimmed
 * @param metrics the Metrics instance to update the metrics
 * @param raw_length
 * @param mgr the FlowManager instance to access the adapter and SID pool information
 * @param trim_info the TrimInfo object containing the results of the demuxing and trimming for the read
 **/
template <typename TrimInfo>
void Demux::RunSimplexAdapter(FixedReadRecord& read, SimplexMetrics& metrics, u32 raw_length, const FlowManager& mgr,
                              const TrimInfo& trim_info) {
  ++metrics.preassign_passing_read_count;
  if (trim_info.insert.Length() < mgr.param.min_trimmed_read_len) {
    // read is too short after trimming, skip writing
    ++metrics.too_short_trimmed_read_count;
    read.SetStatus(FixedReadRecord::Status::kTrimmedTooShortFail);
    return;
  }

  metrics.RecordReadMetrics(trim_info, raw_length);

  if (trim_info.sid.has_value()) {
    // Check if discordant reads should be discarded based on the configured mode
    if (ShouldDiscardDiscordant(trim_info, mgr.param.discordant_sid_mode)) {
      ++metrics.sid_discordant_discarded_read_count;
      read.SetStatus(FixedReadRecord::Status::kDiscordantSidFail);
      return;
    }
    const auto sid{trim_info.sid.value()};
    if (sid >= mgr.sid_pool.size()) [[unlikely]] {
      throw error::Error("SID {} from LUT exceeds SidPool size {}", sid, mgr.sid_pool.size());
    }
    // In kAllSplit mode, partial reads are routed to a separate set of sinks
    // at offset sid_pool.size() from the normal (full) sinks.
    const auto partial_offset = static_cast<u32>(
        (mgr.param.read_length_mode == ReadLengthMode::kAllSplit && IsPartialRead(trim_info)) ? mgr.sid_pool.size()
                                                                                              : 0);
    read.file_sid = 1 + sid + partial_offset;
    read.SetStatus(FixedReadRecord::Status::kDemultiplexed);
  }
}

/**
 * @brief RunDuplexAdapter is a helper function to run the demuxing for duplex adapters.
 * @param read the read to be demultiplexed
 * @param duplex_metrics the DuplexMetrics instance to update the metrics
 * @param mgr the FlowManager instance to access the adapter and SID pool information
 **/
void Demux::RunDuplexAdapter(FixedReadRecord& read, DuplexMetrics& duplex_metrics, const FlowManager& mgr) {
  // actual demuxing done here
  mgr.HairpinFinderObj().FindHairpin(read, duplex_metrics);

  if (read.trim_info_duplex.duplex_status == TrimInfoDuplex::DuplexStatus::kMidAdapterFound) {
    read.file_sid = 1 + read.trim_info_duplex.matches[DuplexMatch::kSID3p].match.barcode_id;
    read.SetStatus(FixedReadRecord::Status::kDemultiplexed);
  } else {
    // all other states fail
    read.SetStatus(FixedReadRecord::Status::kDuplexMidAdapterFail);
  }
}

void Demux::operator()() {
  try {
    StopWatch sw;
    auto& mgr{context.GetManager()};
    auto& batch{context.GetBatchData(batch_nr)};
    auto& metrics{SimplexMetrics::Instance()};
    auto& duplex_metrics{DuplexMetrics::Instance()};
    if (batch.records) {
      for (size_t i = 0; i < batch.Size(); ++i) {
        auto& read = (*(batch.records))[i];
        // filter and record raw reads here
        if (mgr.is_duplex) {
          duplex_metrics.total_length_distr.AddCountToHistogram(read.SeqLen(), 1);
          if (read.SeqLen() < mgr.param.min_read_len) {
            read.SetStatus(FixedReadRecord::Status::kTooShortFail);
            duplex_metrics.unassigned_counts.read_too_short += 1;
            duplex_metrics.unassigned_counts.raw_bases += read.SeqLen();
            duplex_metrics.unassigned_length_distr.AddCountToHistogram(read.SeqLen(), 1);
            continue;
          }
          if (read.GetStatus() == FixedReadRecord::Status::kTooLongFail) {
            duplex_metrics.unassigned_counts.read_too_long += 1;
            duplex_metrics.unassigned_counts.raw_bases += read.SeqLen();
            duplex_metrics.unassigned_length_distr.AddCountToHistogram(read.SeqLen(), 1);
            continue;
          }
        } else {
          metrics.IncrementInputReadCount();
          if (read.SeqLen() < mgr.param.min_read_len) {
            ++metrics.too_short_read_count;
            read.SetStatus(FixedReadRecord::Status::kTooShortFail);
            // read is too short, skip read
            continue;
          }
          if (read.GetStatus() == FixedReadRecord::Status::kTooLongFail) {
            // read is too long, skip read
            continue;
          }
        }

        // Demux and trim the read according to the adapter type
        switch (mgr.adapter_type) {
          using enum AdapterType;
          case kYsu: {
            const auto raw_length = read.SeqLen();
            // actual demuxing done here
            read.trim_info_ysu = mgr.DemuxObjectYsu()(read);
            RunSimplexAdapter(read, metrics, raw_length, mgr, read.trim_info_ysu);
            break;
          }
          case kYs: {
            const auto raw_length = read.SeqLen();
            // actual demuxing done here
            read.trim_info_simplex = mgr.DemuxObjectYs()(read);
            RunSimplexAdapter(read, metrics, raw_length, mgr, read.trim_info_simplex);
            break;
          }
          case kSimplex: {
            const auto raw_length = read.SeqLen();
            RunSimplexAdapter(read, metrics, raw_length, mgr, mgr.DemuxObjectSimplex()(read));
            break;
          }
          case kSimplex10x: {
            const auto raw_length = read.SeqLen();
            read.trim_info_simplex = mgr.DemuxObjectSimplex10x()(read);
            RunSimplexAdapter(read, metrics, raw_length, mgr, read.trim_info_simplex);
            break;
          }
          case kDuplex:
          case kDuplexUMI:
          case kDuplexStem:
            RunDuplexAdapter(read, duplex_metrics, mgr);
            break;
          default:
            throw error::Error("Unknown adapter type found. Saw adapter type: {}", static_cast<u32>(mgr.adapter_type));
        }
      }
    }

    context.AddToDemuxTime(sw.ElapsedTime());
  } catch (const std::exception& e) {
    Logging::Error("Demux::operator() failed: {}", e.what());
    SetTaskException(std::current_exception());
  }
}
}  // namespace xoos::demux
