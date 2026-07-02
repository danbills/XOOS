#pragma once
#include <xoos/error/error.h>

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <limits>
#include <span>
#include <string_view>
#include <system_error>

#include "io/read-record.h"
#include "task/task.h"
#include "utility/alignment-util.h"

namespace xoos::demux {

constexpr char kOriginalIdPrefix = '@';
constexpr char kReadNameSeparator = ':';
constexpr char kMissingBarcodeChar = '*';
constexpr char kQualitySeparator = '+';
constexpr char kBitFlagSidSeparator = '|';
constexpr char kCommentSeparator = ' ';
constexpr std::string_view kTrimDebugRawLen = " rl:i:";
constexpr std::string_view kTrimDebugTrim5p = " ta:i:";
constexpr std::string_view kTrimDebugTrim3p = " tb:i:";

// FASTQ header bitflag layout: a single byte encoded as a decimal integer between the read name and
// the SID, e.g. `@<name>:<bitflag>|<sid>`. A set bit means the corresponding feature is present.
// 5' SID-present bit (0b1000) in the FASTQ header bitflag.
inline constexpr auto kBitFlagSid5p = static_cast<std::byte>(0b1000);
// 3' SID-present bit (0b0100) in the FASTQ header bitflag.
inline constexpr auto kBitFlagSid3p = static_cast<std::byte>(0b0100);
// 5' UMI-present bit (0b0010) in the FASTQ header bitflag.
inline constexpr auto kBitFlagUmi5p = static_cast<std::byte>(0b0010);
// 3' UMI-present bit (0b0001) in the FASTQ header bitflag.
inline constexpr auto kBitFlagUmi3p = static_cast<std::byte>(0b0001);
// Convenience mask: duplex reads always carry both SIDs.
inline constexpr auto kBitFlagDuplexSids = kBitFlagSid5p | kBitFlagSid3p;

// Maximum number of decimal digits required to render a u8 (the FASTQ header bitflag).
inline constexpr size_t kMaxU8DecimalDigits = std::numeric_limits<u8>::digits10 + 1;
// Maximum number of decimal digits required to render a u32 (SID, raw length, trim positions).
inline constexpr size_t kMaxU32DecimalDigits = std::numeric_limits<u32>::digits10 + 1;

/**
 * @brief Options controlling how simplex reads are formatted into FASTQ output.
 *
 * Bundled together so we don't have to thread multiple positional bools through the formatting call chain.
 */
struct SimplexFormatOptions {
  /// If true, preserve original quality scores instead of overriding with kSimplexBaseQual.
  bool suppress_qual_override = false;
  /// If true, append trim coordinates (rl/ta/tb tags) to the FASTQ comment.
  bool trimming_debug = false;
};

/**
 * @brief Writes demultiplexed data into appropriate output files.
 *
 * The Formatter class converts read records into various output formats for each adapter and writes them to disk.
 * It operates as a task within the taskflow scheduler framework.
 */
class Formatter : public Task {
 public:
  /**
   * @brief Constructs a Formatter task.
   * @param[in] exec The flow context providing execution environment.
   * @param[in] batch_nr The batch number to process.
   */
  Formatter(FlowContext& exec, size_t batch_nr);
  ~Formatter() override = default;

  /**
   * @brief Executes the formatter task.
   *
   * Called by the task scheduler to write a batch of records from the demux data into the appropriate output files.
   */
  void operator()();

  /**
   * @brief Converts a FixedReadRecord to processed FASTQ format for duplex adapter data. Outputs processed format
   * including YC tag, consensus sequence, etc.
   *
   * Exposed as static for projects that use demux like an API in another context.
   *
   * @param read The read record to convert.
   * @param p_data Output buffer for the formatted data.
   * @param offset Current position in buffer; updated after write.
   * @details Assumes the FixedReadRecord provides any additional metadata (for example a YC tag) via its
   * comment/header fields, and that ConsensusSeq() contains the duplex consensus sequence.
   */
  static void ToFastqDuplex(const FixedReadRecord& read, char* p_data, size_t& offset);

  /**
   * @brief Converts a FixedReadRecord to original raw FASTQ format. Outputs the original unprocessed read data. Usually
   * used for failed reads.
   *
   * Exposed as static for projects that use demux like an API in another context.
   *
   * @param[in] read The read record to convert.
   * @param[out] p_data Output buffer for the formatted data.
   * @param[in,out] offset Current position in buffer; updated after write.
   * @warning Assumes Seq(), Qual(), Name(), Comment() are unaltered, which may not hold true due to partial processing
   */
  static void ToFastqRaw(const FixedReadRecord& read, char* p_data, size_t& offset);

 private:
  /**
   * @brief Determine whether a record should be skipped during the initial scan for unwritten records.
   *
   * @param record The record to evaluate.
   * @return true if the record should be skipped (already processed or not eligible).
   */
  bool ShouldSkipRecord(const FixedReadRecord& record) const;

  /**
   * @brief Format a single record into the output buffer, dispatching on its status.
   *
   * Encapsulates the status-based switch logic previously inlined in operator().
   *
   * @param[in,out] record    The record to format; status is updated to kWritten on success.
   * @param         sid       The SID currently being collected.
   * @param[in,out] mem       The output buffer.
   * @return true if bytes were written to mem (caller should bounds-check), false if skipped.
   */
  bool FormatRecord(FixedReadRecord& record, u32 sid, FormattedOutput& mem) const;

  /**
   * @brief Schedule async write tasks for the collected sinks and mark the batch complete.
   *
   * @param sinks Span of sink descriptors populated during formatting.
   */
  void ScheduleWriteTasks(std::span<SinkData> sinks) const;

  /**
   * @brief Write a successfully-demultiplexed record to the buffer, dispatching on the configured adapter type.
   *
   * Reads the simplex formatting options and adapter type from the enclosing flow context, so they do not have to
   * be threaded through from `operator()`.
   */
  void PassingReadToBuffer(FixedReadRecord& read, FormattedOutput& mem) const;

  /**
   * @brief Write a failed record to the buffer as raw FASTQ.
   *
   * @warning Assumes data is unaltered, which may not hold true due to partial processing.
   */
  static void FailedReadToBuffer(FixedReadRecord& read, FormattedOutput& mem);
};

void WriteUintToBuffer(char* p_data, size_t& offset, size_t buffer_end, u32 value, std::string_view error_msg);
void WriteUmiValue(char* p_data, size_t& offset, size_t buffer_end, const std::optional<u32>& umi_value);
void WriteFastqHeaderWithBitflagAndSid(const FixedReadRecord& read, char* p_data, size_t& offset, std::byte bitflag);

/**
 * @brief Append a contiguous byte range to the output buffer and advance the offset.
 *
 * For trivially-copyable contiguous ranges the standard library lowers `std::copy_n` to
 * `__builtin_memmove`, which on modern libc shares fast paths with `memcpy`; functionally
 * identical here because the destination never overlaps the source. Caller is responsible
 * for ensuring the destination has sufficient capacity.
 */
inline void AppendBytes(char* const p_data, size_t& offset, const std::string_view src) {
  std::copy_n(src.data(), src.size(), p_data + offset);
  offset += src.size();
}

/// @brief Append trim debug tags to the FASTQ comment: " rl:i:<raw_len> ta:i:<spos> tb:i:<epos>"
void AppendTrimDebugTags(char* p_data, size_t& offset, size_t buffer_end, u32 raw_len, const LociRange& insert);

void ToFastqDuplexUMI(const FixedReadRecord& read, FormattedOutput& mem);

/**
 * @brief Write a simplex read to the buffer without UMI information in the read name.
 *
 * @tparam TrimInfo The type containing trimming information (e.g., TrimInfoYsu, TrimInfoSimplex).
 * @param read The read to write.
 * @param mem The buffer to write to.
 * @param trim_info The trimming information for the read, used to extract the trimmed sequence if applicable.
 * @param opts Formatting options (quality override, trim debug tags).
 *
 * @details Output format: @<original id> <original comment>
 *          This function writes the read name without any UMI sequences, followed by the comment,
 *          sequence (trimmed if applicable), quality separator, and quality scores.
 *          Used for adapter types that don't capture UMI information (e.g., kYs).
 */
template <typename TrimInfo>
void ToFastqSimplex(const FixedReadRecord& read, FormattedOutput& mem, const TrimInfo& trim_info,
                    SimplexFormatOptions opts) {
  auto& offset = mem.nr_bytes;
  auto* const p_data = mem.p_data;

  // Write header: @<name>:<bitflag>|<sid>
  auto bitflag = GenerateSidBitFlagFromTrimInfo(trim_info);
  WriteFastqHeaderWithBitflagAndSid(read, p_data, offset, bitflag);

  // Write the body (comment, sequence, quality scores, etc.)
  AppendFastqSimplexBody(read, p_data, offset, trim_info, opts, mem.capacity);
}

/**
 * @brief Writes the read body (comment, sequence, and quality scores) to the buffer,
 *        assuming the read name has already been written.
 * Useful to reduce code duplication between simplex UMI and non-UMI output.
 *
 * @param read The read to write.
 * @param p_data Output buffer for the formatted data.
 * @param offset Current position in buffer; updated after write.
 * @param trim_info The trimming and UMI information for the read.
 * @param opts Formatting options (quality override, trim debug tags).
 * @param buffer_end End of the output buffer (for bounds checking in trim debug tags).
 */
template <typename TrimInfo>
void AppendFastqSimplexBody(const FixedReadRecord& read, char* const p_data, size_t& offset, const TrimInfo& trim_info,
                            const SimplexFormatOptions opts, const size_t buffer_end) {
  // Add space between comment and name
  p_data[offset++] = ' ';

  // Copy the comment as-is - and don't forget to add the '\n' at the end.
  AppendBytes(p_data, offset, {read.Comment(), read.CommentLen()});

  // Append trim debug tags: rl:i:<raw_len> ta:i:<insert_spos> tb:i:<insert_epos>
  if (opts.trimming_debug) {
    AppendTrimDebugTags(p_data, offset, buffer_end, read.SeqLen(), trim_info.insert);
  }

  // Add a newline after the comment.
  p_data[offset++] = '\n';

  // Determine the source pointers and length based on whether the read was trimmed.
  const auto insert_offset = trim_info.sid ? trim_info.insert.spos : 0;
  const auto length = trim_info.sid ? trim_info.insert.Length() : read.SeqLen();

  // Copy the sequence.
  AppendBytes(p_data, offset, {read.Seq() + insert_offset, length});

  p_data[offset++] = '\n';
  p_data[offset++] = kQualitySeparator;
  p_data[offset++] = '\n';

  // Write quality scores. By default, override with kSimplexBaseQual.
  // When suppress_qual_override is true, preserve the original quality scores.
  if (opts.suppress_qual_override) {
    AppendBytes(p_data, offset, {read.Qual() + insert_offset, length});
  } else {
    std::ranges::fill_n(p_data + offset, static_cast<std::ptrdiff_t>(length), kSimplexBaseQual);
    offset += length;
  }

  // Add a newline at the end of the quality scores.
  p_data[offset++] = '\n';
}

/**
 * @brief Write the read to the buffer, using the TrimInfo to get UMI sequences and trimming information.
 * @param read The read to write.
 * @param mem The buffer to write to.
 * @param trim_info The trimming and UMI information for the read.
 * @param opts Formatting options (quality override, trim debug tags).
 */
template <typename TrimInfo>
void ToFastqSimplexUMI(const FixedReadRecord& read, FormattedOutput& mem, const TrimInfo& trim_info,
                       SimplexFormatOptions opts) {
  auto& offset = mem.nr_bytes;
  auto* const p_data = mem.p_data;

  // Write header: @<name>:<bitflag>|<sid>
  auto bitflag = GenerateUmiBitFlagFromTrimInfo(trim_info);
  WriteFastqHeaderWithBitflagAndSid(read, p_data, offset, bitflag);

  // UMI name portion, parsers look for a pipe '|' character currently
  p_data[offset++] = kReadNameSeparator;
  if (trim_info.umi_5p.has_value()) {
    auto [ptr, ec] = std::to_chars(p_data + offset, p_data + mem.capacity, *trim_info.umi_5p);
    offset = static_cast<size_t>(ptr - p_data);
    if (ec != std::errc()) {
      throw error::Error("UMI conversion error. Either buffer too small or value invalid.\n{}",
                         std::make_error_code(ec).message());
    }
  } else {
    p_data[offset++] = kMissingBarcodeChar;
  }
  p_data[offset++] = kReadNameSeparator;
  if (trim_info.umi_3p.has_value()) {
    auto [ptr, ec] = std::to_chars(p_data + offset, p_data + mem.capacity, *trim_info.umi_3p);
    offset = static_cast<size_t>(ptr - p_data);
    if (ec != std::errc()) {
      throw error::Error("UMI conversion error. Either buffer too small or value invalid.\n{}",
                         std::make_error_code(ec).message());
    }
  } else {
    p_data[offset++] = kMissingBarcodeChar;
  }

  // Write the body (comment, sequence, quality scores, etc.)
  AppendFastqSimplexBody(read, p_data, offset, trim_info, opts, mem.capacity);
}

/**
 * @brief Generates a SID-only bitflag indicating which SIDs are present.
 * UMI bits are always zero; use GenerateUmiBitFlagFromTrimInfo when UMI presence matters.
 * @tparam TrimInfo Any trim-info type with optional sid_5p / sid_3p fields
 *         (e.g., TrimInfoYsu, TrimInfoSimplex).
 * @param trim_info The trimming information to check.
 * @return Bitflag: 5' sid (0b1000), 3' sid (0b0100), 5' umi (0b0010), 3' umi (0b0001)
 */
template <typename TrimInfo>
std::byte GenerateSidBitFlagFromTrimInfo(const TrimInfo& trim_info) {
  auto bitflag = std::byte{0};
  if (trim_info.sid_5p.has_value()) {
    bitflag |= kBitFlagSid5p;
  }
  if (trim_info.sid_3p.has_value()) {
    bitflag |= kBitFlagSid3p;
  }
  return bitflag;
}

/**
 * @brief Generates a bitflag based on which SIDs and UMIs are present in the given TrimInfo.
 * @tparam TrimInfo Any trim-info type with optional sid_5p / sid_3p / umi_5p / umi_3p fields
 *         (e.g., TrimInfoYsu, TrimInfoDuplexUMI).
 * @param trim_info The trimming and UMI information to check.
 * @return Bitflag: 5' sid (0b1000), 3' sid (0b0100), 5' umi (0b0010), 3' umi (0b0001)
 */
template <typename TrimInfo>
std::byte GenerateUmiBitFlagFromTrimInfo(const TrimInfo& trim_info) {
  auto bitflag = std::byte{0};
  if (trim_info.sid_5p.has_value()) {
    bitflag |= kBitFlagSid5p;
  }
  if (trim_info.sid_3p.has_value()) {
    bitflag |= kBitFlagSid3p;
  }
  if (trim_info.umi_5p.has_value()) {
    bitflag |= kBitFlagUmi5p;
  }
  if (trim_info.umi_3p.has_value()) {
    bitflag |= kBitFlagUmi3p;
  }
  return bitflag;
}

std::byte GenerateUmiBitFlagFromFixedDuplexRead(const FixedReadRecord& read);

}  // namespace xoos::demux
