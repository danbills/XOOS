#include "task/formatter.h"

#include <fmt/format.h>
#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include <algorithm>
#include <array>
#include <span>

#include "task/flow-context.h"
#include "task/flow-manager.h"
#include "utility/stop-watch.h"

namespace xoos::demux {
Formatter::Formatter(FlowContext& exec, size_t batch_nr) : Task(exec, batch_nr, fmt::format("Writer{}", batch_nr)) {}

/**
 * @brief Generates a bitflag for the read based on which UMIs are present in the TrimInfoDuplexUMI. The bitflag is used
 * to indicate the presence of 5' and 3' SIDs and UMIs in the output FASTQ header for duplex reads.
 * @param[in] read The read to generate the bitflag for. We will look at the TrimInfoDuplexUMI to determine which UMIs
 * were present and set the appropriate bits in the bitflag.
 * @return
 */
std::byte GenerateUmiBitFlagFromFixedDuplexRead(const FixedReadRecord& read) {
  auto bitflag = kBitFlagDuplexSids;
  if (read.trim_info_duplex.umi_5p.has_value()) {
    bitflag |= kBitFlagUmi5p;
  }
  if (read.trim_info_duplex.umi_3p.has_value()) {
    bitflag |= kBitFlagUmi3p;
  }
  return bitflag;
}

/**
 * @brief Writes the read body (YC tag, comment, consensus sequence, and quality scores) to the buffer,
 *        assuming the read name has already been written.
 * Useful to reduce code duplication between UMI and non-UMI output.
 *
 * @param read The read to write.
 * @param p_data Output buffer for the formatted data.
 * @param offset Current position in buffer; updated after write.
 */
static void AppendFastqDuplexBody(const FixedReadRecord& read, char* const p_data, size_t& offset) {
  // Add the YC tag
  AppendBytes(p_data, offset, {read.Seq(), read.iupac_length});

  // Copy the comment as-is - and don't forget to add the '\n' at the end.
  AppendBytes(p_data, offset, {read.Comment(), read.CommentLen()});
  // Add a newline after the comment.
  p_data[offset++] = '\n';

  const auto consensus_len = static_cast<size_t>(read.consensus_seq_len);

  // The sequence is now a consensus sequence.
  AppendBytes(p_data, offset, {read.ConsensusSeq(), consensus_len});

  // Add a newline at the end of the sequence.
  p_data[offset++] = '\n';
  // Add the '+' character after the sequence.
  p_data[offset++] = kQualitySeparator;
  // Add a newline after the '+' character.
  p_data[offset++] = '\n';

  // Copy the quality scores; these artificial scores have replaced the sequence.
  AppendBytes(p_data, offset, {read.Qual(), consensus_len});

  // Add a newline at the end of the quality scores.
  p_data[offset++] = '\n';
}

/**
 * @brief Write the read to the buffer, using the DuplexUmi bundle to get the UMI sequences.
 * @param read The read to write.
 * @param mem The buffer to write to.
 * @details Implementation is similar to the Duplex version, but we add the UMI sequences to the read name.
 * '*' is used to indicate missing barcodes.
 */
void ToFastqDuplexUMI(const FixedReadRecord& read, FormattedOutput& mem) {
  auto& offset = mem.nr_bytes;
  auto* const p_data = mem.p_data;

  // Write header: @<name>:<bitflag>|<sid>
  const auto bitflag = GenerateUmiBitFlagFromFixedDuplexRead(read);
  WriteFastqHeaderWithBitflagAndSid(read, p_data, offset, bitflag);

  // Write UMI values: :<umi_5p>:<umi_3p>
  p_data[offset++] = kReadNameSeparator;
  WriteUmiValue(p_data, offset, mem.capacity, read.trim_info_duplex.umi_5p);
  p_data[offset++] = kReadNameSeparator;
  WriteUmiValue(p_data, offset, mem.capacity, read.trim_info_duplex.umi_3p);

  AppendFastqDuplexBody(read, p_data, offset);
}

/**
 * @brief Write original read to buffer. Intended for failed reads.
 * @param read The read to write.
 * @param mem The buffer to write to.
 * @warning This function assumes data is unaltered, which may not hold true due to partial processing.
 */
void Formatter::FailedReadToBuffer(FixedReadRecord& read, FormattedOutput& mem) {
  auto& offset = mem.nr_bytes;
  auto* const p_data = mem.p_data;

  Formatter::ToFastqRaw(read, p_data, offset);
  read.SetStatus(FixedReadRecord::Status::kWritten);
}

/**
 * @brief Helper function to write a uint value to buffer using std::to_chars with error handling.
 * @param p_data Output buffer for the formatted data.
 * @param offset Current position in buffer; updated after write.
 * @param buffer_end End of the buffer.
 * @param value The value to write.
 * @param error_msg Error message prefix if conversion fails.
 */
void WriteUintToBuffer(char* const p_data, size_t& offset, const size_t buffer_end, const u32 value,
                       const std::string_view error_msg) {
  const auto [ptr, ec] = std::to_chars(p_data + offset, p_data + buffer_end, value);
  if (ec != std::errc()) {
    throw error::Error("{} Either buffer too small or value invalid.\n", error_msg, std::make_error_code(ec).message());
  }
  offset = static_cast<size_t>(ptr - p_data);
}

/**
 * @brief Helper function to write an optional UMI value to buffer.
 * Writes the UMI value if present, otherwise writes '*'.
 * @param p_data Output buffer for the formatted data.
 * @param offset Current position in buffer; updated after write.
 * @param buffer_end End of the buffer.
 * @param umi_value Optional UMI value to write.
 */
void WriteUmiValue(char* const p_data, size_t& offset, const size_t buffer_end, const std::optional<u32>& umi_value) {
  if (umi_value.has_value()) {
    WriteUintToBuffer(p_data, offset, buffer_end, *umi_value, "UMI conversion error.");
  } else {
    p_data[offset++] = kMissingBarcodeChar;
  }
}

/**
 * @brief Helper function to write the FASTQ header with bitflag and SID.
 * Writes: @<name>:<bitflag>|<sid>
 * @param read The read record.
 * @param p_data Output buffer for the formatted data.
 * @param offset Current position in buffer; updated after write.
 * @param bitflag The bitflag value to write.
 */
void WriteFastqHeaderWithBitflagAndSid(const FixedReadRecord& read, char* const p_data, size_t& offset,
                                       const std::byte bitflag) {
  // Write '@' and name
  p_data[offset++] = kOriginalIdPrefix;
  AppendBytes(p_data, offset, {read.Name(), read.NameLen()});

  // Write ':' and bitflag
  p_data[offset++] = kReadNameSeparator;
  WriteUintToBuffer(p_data, offset, offset + kMaxU8DecimalDigits, static_cast<u32>(bitflag),
                    "Bitflag conversion error.");

  // Write '|' and SID.
  p_data[offset++] = kBitFlagSidSeparator;
  WriteUintToBuffer(p_data, offset, offset + kMaxU32DecimalDigits, read.file_sid, "SID conversion error.");
}

/**
 * @brief Append trim debug tags to the FASTQ comment.
 * Format: " rl:i:<raw_len> ta:i:<spos> tb:i:<epos>"
 * @param p_data Output buffer.
 * @param offset Current position in buffer; updated after write.
 * @param buffer_end End of the output buffer for bounds checking.
 * @param raw_len The raw (untrimmed) read length.
 * @param insert The insert range from trimming.
 */
void AppendTrimDebugTags(char* const p_data, size_t& offset, const size_t buffer_end, const u32 raw_len,
                         const LociRange& insert) {
  // Worst case: 3 tag prefixes + 3 u32 values rendered in decimal.
  constexpr size_t kMaxTrimDebugBytes =
      kTrimDebugRawLen.size() + kTrimDebugTrim5p.size() + kTrimDebugTrim3p.size() + 3 * kMaxU32DecimalDigits;
  if (offset + kMaxTrimDebugBytes > buffer_end) {
    throw error::Error("Trim debug tags require {} bytes but only {} available.", kMaxTrimDebugBytes,
                       buffer_end - offset);
  }

  // Write rl:i:<raw_len>
  AppendBytes(p_data, offset, kTrimDebugRawLen);
  WriteUintToBuffer(p_data, offset, buffer_end, raw_len, "Failed to write trim debug raw length.");

  // Write ta:i:<spos>
  AppendBytes(p_data, offset, kTrimDebugTrim5p);
  WriteUintToBuffer(p_data, offset, buffer_end, insert.spos, "Failed to write trim debug 5' position.");

  // Write tb:i:<epos>
  AppendBytes(p_data, offset, kTrimDebugTrim3p);
  WriteUintToBuffer(p_data, offset, buffer_end, insert.epos, "Failed to write trim debug 3' position.");
}

/**
 * @brief Write the record to the buffer based on the configured adapter type.
 * Reads adapter type and simplex formatting options from the enclosing flow context.
 * @param read The read record to write.
 * @param mem The buffer to write to.
 */
void Formatter::PassingReadToBuffer(FixedReadRecord& read, FormattedOutput& mem) const {
  using enum AdapterType;
  const auto& params = context.DemuxParam();
  const SimplexFormatOptions simplex_opts{.suppress_qual_override = params.suppress_simplex_qual_override,
                                          .trimming_debug = params.trimming_debug};
  switch (context.GetManager().adapter_type) {
    case kYsu:
      ToFastqSimplexUMI(read, mem, read.trim_info_ysu, simplex_opts);
      break;
    case kYs:
    case kSimplex:
    case kSimplex10x:
      // These use the same TrimInfoSimplex layout and FASTQ formatting logic.
      ToFastqSimplex(read, mem, read.trim_info_simplex, simplex_opts);
      break;
    case kDuplexStem:
    case kDuplex:
      Formatter::ToFastqDuplex(read, mem.p_data, mem.nr_bytes);
      if (mem.nr_bytes > mem.capacity) {
        throw error::Error("Sequence too long to store in internal buffer!");
      }
      break;
    case kDuplexUMI:
      ToFastqDuplexUMI(read, mem);
      break;
    default:
      throw error::Error("Unknown adapter type {} encountered in Formatter",
                         static_cast<s32>(context.GetManager().adapter_type));
  }
  read.SetStatus(FixedReadRecord::Status::kWritten);
}

void Formatter::ToFastqDuplex(const FixedReadRecord& read, char* const p_data, size_t& offset) {
  // Write header: @<name>:<bitflag>|<sid>
  const auto bitflag = GenerateUmiBitFlagFromFixedDuplexRead(read);
  WriteFastqHeaderWithBitflagAndSid(read, p_data, offset, bitflag);

  // Write the body (YC tag, comment, consensus sequence, quality scores, etc.)
  AppendFastqDuplexBody(read, p_data, offset);
}

void Formatter::ToFastqRaw(const FixedReadRecord& read, char* const p_data, size_t& offset) {
  p_data[offset++] = kOriginalIdPrefix;
  AppendBytes(p_data, offset, {read.Name(), read.NameLen()});
  // copy the comment
  p_data[offset++] = kCommentSeparator;
  AppendBytes(p_data, offset, {read.Comment(), read.CommentLen()});
  // Add a newline after the name.
  p_data[offset++] = '\n';

  // Copy the sequence.
  AppendBytes(p_data, offset, {read.Seq(), read.SeqLen()});

  // Separator line between the sequence and the quality scores.
  p_data[offset++] = '\n';
  p_data[offset++] = kQualitySeparator;
  p_data[offset++] = '\n';

  // Copy the quality scores.
  AppendBytes(p_data, offset, {read.Qual(), read.QualLen()});
  // Add a newline after the quality scores
  p_data[offset++] = '\n';
}

bool Formatter::ShouldSkipRecord(const FixedReadRecord& record) const {
  if (context.DemuxParam().output_failed_reads) {
    return record.GetStatus() == FixedReadRecord::Status::kWritten;
  }
  return record.GetStatus() != FixedReadRecord::Status::kDemultiplexed;
}

bool Formatter::FormatRecord(FixedReadRecord& record, const u32 sid, FormattedOutput& mem) const {
  auto& mgr{context.GetManager()};
  using enum FixedReadRecord::Status;

  switch (record.GetStatus()) {
    case kNotRead:
      throw error::Error("Internal error: record in Formatter with status {} encountered.",
                         static_cast<s32>(record.GetStatus()));
    case kWritten:
      // Nothing to do here, but skip buffer checking
    case kTooLongFail:
      // The read being too long means we don't have it in the buffer, so we can't write it out.
      // TODO: Find some way of preserving this read
      return false;
    case kRead:
      // shouldn't be possible but we don't always track status for simplex reads and will be considered failed
      // TODO: Make simplex adapters use failure states
    case kTooShortFail:
    case kTrimmedTooShortFail:
    case kDuplexMidAdapterFail:
    case kDuplexEditDistanceFail:
    case kDuplexTooLongFail:
    case kFailedMidadapterTrimFail:
    case kDiscordantSidFail:
      // Failed reads are only written when --output-failed-reads is set and the read is unassigned.
      if (!context.DemuxParam().output_failed_reads || sid != FixedReadRecord::kUnassignedSID) {
        return false;
      }
      // For simplex adapters, override quality scores in-place before writing raw output.
      if (!mgr.is_duplex && !context.DemuxParam().suppress_simplex_qual_override) {
        std::ranges::fill_n(record.Qual(), record.QualLen(), kSimplexBaseQual);
      }
      FailedReadToBuffer(record, mem);
      return true;
    case kDemultiplexed:
      if (sid != record.file_sid) {
        return false;
      }
      PassingReadToBuffer(record, mem);
      return true;
    default:
      throw error::Error("Unknown FixedReadRecord status {} encountered in formatter.",
                         static_cast<s32>(record.GetStatus()));
  }
}

void Formatter::ScheduleWriteTasks(const std::span<SinkData> sinks) const {
  auto& executor{context.GetManager().Executor()};

  if (!sinks.empty()) {
    if (sinks.size() > FormattedOutput::kMaxNumberSinks) {
      throw error::Error("Too many sinks. Expected at most {} sinks, but got {} sinks.",
                         FormattedOutput::kMaxNumberSinks, sinks.size());
    }

    // Now create the write tasks. These tasks will not need to be stored, so creating them on the stack.
    std::array<tf::AsyncTask, FormattedOutput::kMaxNumberSinks> write_tasks{};
    for (size_t i = 0; i < sinks.size(); ++i) {
      write_tasks[i] = context.CreateWriteTask(sinks[i]);
    }
    context.nr_writes_scheduled += sinks.size();

    // Join task: completes once all per-sink write_tasks finish; no work of its own.
    context.SetWriteTask(batch_nr,
                         executor.silent_dependent_async([] { /* barrier task, no work */ }, write_tasks.begin(),
                                                         write_tasks.begin() + sinks.size()));
  } else {
    // No sinks created, so we need to create a dummy write task.
    context.SetWriteTask(batch_nr, executor.silent_dependent_async([] { /* no-op placeholder write task */ }));
  }
}

void Formatter::operator()() {
  try {
    const StopWatch sw;
    auto& batch{context.GetBatchData(batch_nr)};

    size_t nr_sinks = 0;
    size_t current_index = 0;
    auto& mem{batch.formatted_output};
    mem.nr_bytes = 0;
    auto* const sink_data = batch.formatted_output.file_sinks;

    while (current_index != batch.num_records) {
      // Advance past records that are already handled or not eligible.
      while (current_index < batch.num_records && ShouldSkipRecord((*batch.records)[current_index])) {
        ++current_index;
      }
      if (current_index >= batch.num_records) {
        break;
      }

      if (nr_sinks >= FormattedOutput::kMaxNumberSinks) {
        throw error::Error(
            "Error condition detected: max number SIDs in a batch ({}) exceeded, generated output would be invalid.\n"
            "Saw {} sinks so far.",
            FormattedOutput::kMaxNumberSinks, nr_sinks);
      }

      const auto sid = (*batch.records)[current_index].file_sid;
      const auto start_offset = mem.nr_bytes;
      auto& this_sink = sink_data[nr_sinks];
      this_sink.p_data = mem.p_data + start_offset;
      this_sink.sid_id = sid;
      this_sink.batch_id = batch_nr;
      this_sink.flow_context = &context;

      // Format all records matching this SID into the buffer.
      for (auto index = current_index; index < batch.num_records; ++index) {
        if (!FormatRecord((*batch.records)[index], sid, mem)) {
          continue;
        }
        // Written something — bounds-check the buffer.
        if (start_offset + mem.nr_bytes > mem.capacity) {
          throw error::Error("Write buffer not large enough. Bytes available: {}, bytes needed: {}", mem.capacity,
                             start_offset + mem.nr_bytes);
        }
      }

      this_sink.length = mem.nr_bytes - start_offset;
      ++nr_sinks;
    }

    ScheduleWriteTasks({sink_data, nr_sinks});

    // Mark batch as completed, but this also schedules a new batch if input is available.
    context.MarkBatchCompleted();
    context.AddToFormatTime(static_cast<u64>(sw.ElapsedTime()));
  } catch (const std::exception& e) {
    Logging::Error("Formatter::operator() failed: {}", e.what());
    SetTaskException(std::current_exception());
  }
}
}  // namespace xoos::demux
