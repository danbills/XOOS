#include "rdb-sequence-reader.h"

#include <fmt/format.h>
#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/util/parse-int.h>
#include <xoos/util/string-functions.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#include <tuple>

#include "rdb-2bit-utils.h"
#include "utility/math-util.h"
#include "xoos/enum/enum-util.h"

namespace xoos::demux {

// Constants for RDB file format
constexpr u32 kAppVersionOffset = 8;
constexpr u32 kAppVersionLength = 32;
constexpr u32 kRunIdOffset = 40;
constexpr u32 kRunIdLength = 64;
constexpr u32 kBatchIndexOffset = 104;
constexpr u32 kBlockTypeOffset = 110;
constexpr u32 kNumDatasetsOffset = 112;
constexpr u32 kTotalDatasetBytesOffset = 128;
constexpr u32 kDatasetSizeOffset = 66;

constexpr u32 kNormalBlock = 0;
constexpr u32 kChunkSize = 4096;
constexpr u32 kBlockInfoSize = 160;
constexpr u32 kDatasetInfoSize = 82;
// Maximum number of datasets that fit in a single RDB block header.
// The header is kChunkSize (4096) bytes; the first kBlockInfoSize (160) bytes are block info,
// and each dataset descriptor is kDatasetInfoSize (82) bytes, giving floor((4096-160)/82) = 48.
// Reading beyond this count would walk past the end of the header buffer.
constexpr u32 kMaxNumDatasets = (kChunkSize - kBlockInfoSize) / kDatasetInfoSize;  // 48 per RDB spec

constexpr u32 kMaxBase64UIDLength = 6;

// Indices into kRequiredDatasetNames for each required dataset.
constexpr size_t kDatasetIdxTracesCount = 0;
constexpr size_t kDatasetIdxTraces = 1;
constexpr size_t kDatasetIdxQscores = 2;
constexpr size_t kDatasetIdxReadLength = 3;
constexpr size_t kDatasetIdxCellIndex = 4;
constexpr size_t kNumRequiredDatasets = 5;

// The 5 datasets that must be present in every normal RDB block.
constexpr std::array<std::string_view, kNumRequiredDatasets> kRequiredDatasetNames = {
    "molecular_traces_count", "molecular_traces", "molecular_traces_qscores", "molecular_traces_read_length",
    "mt_cell_index"};

// Extract a null-terminated string from a fixed-size region of the header buffer.
static std::string ExtractHeaderString(const std::string_view header, const u32 offset, const u32 max_length) {
  if (offset >= header.size()) {
    XOOS_LOG_WARN_ONCE("Offset {} is out of bounds for header of size {}, returning empty header", offset,
                       header.size());
    return {};
  }
  const auto region = header.substr(offset, max_length);
  const auto end_pos = region.find('\0');
  return std::string{end_pos == std::string_view::npos ? region : region.substr(0, end_pos)};
}

// Datapath 2.0 format has 5 sections: date_group_num_instr_cd
// NS3_Gamma1 format has 5 sections.
constexpr u32 kNumPartsDp2OrGamma1 = 5;
// HTP 1.0 format has 6 sections: date_group_num_instr_cd_cycle
// NS3_Gamma2 format has 6 sections. NOTE: order name can have underscores.
constexpr u32 kNumPartsHtp1OrGamma2 = 6;

// if values are not specified by in the prefix queue num
constexpr auto kUnspecifiedQueueNumber = "Q*";

// constant to help logging format
constexpr std::string_view kRunIdTypeSuffix = "-run-id-type";

template <class T>
T GetField(const char* base, u16 offset) {
  return *(reinterpret_cast<const T*>(base + offset));
}

// Validate that all required datasets were found; throw with a detailed message if any are missing.
static void ValidateRequiredDatasets(const std::array<bool, kNumRequiredDatasets>& dataset_found,
                                     const fs::path& file_path, const std::streamoff block_start,
                                     const std::string& rdb_type, const std::string_view header) {
  if (std::ranges::all_of(dataset_found, [](const bool v) { return v; })) {
    return;
  }
  const auto app_version = ExtractHeaderString(header, kAppVersionOffset, kAppVersionLength);
  std::string missing;
  std::string present;
  for (size_t i = 0; i < kRequiredDatasetNames.size(); ++i) {
    auto& target = dataset_found[i] ? present : missing;
    if (!target.empty()) {
      target += ", ";
    }
    target += kRequiredDatasetNames[i];
  }
  throw error::Error(
      "Required datasets missing in '{}' at block offset {}. "
      "RDB type: '{}', app version: '{}'. "
      "Missing: [{}]. Present: [{}]",
      file_path.string(), block_start, rdb_type, app_version, missing, present);
}

/**
 * Here we build the run prefix portion of the read name based on the run_id and cycle_id.
 *
 * The returned prefix will be:
 * {date}:{sequencer}:Q{queue_num}:R{run_num}:{cycle_id}
 *
 *
 *    HTP 1.0 format.
 *    220421_ENG-SYS-HTP_01_cloverleaf_WWX06R05C05_cycle01
 *    0      1           2  3          4            5
 *    date   group       │  instrument chip         cycle
 *                       run
 *   Datapath 2.0 format.
 *    220421_ENG-SYS-HTP_01_cloverleaf_WWX06R05C05
 *    0      1           2  3          4
 *    date   group       │  instrument chip
 *                       run
 *  NS3_Gamma1 format.
 *    20230526_2-00063_13_123456A_1
 *    0        1       2  3       4
 *    date     |       |  tube    cycle
 *             |       run
 *             instrument
 *  NS3_Gamma2 format.
 *    20230526_2-00063_Q13_R12_123456A_1
 *    0        1       2   3   4       5
 *    date     |       |   run order   tube
 *             |       queue
 *             instrument
 *
 * @brief Parse the run_id and determine the RDB run id type.
 *
 * If the run_id doesn't match expected formats, returns the placeholder
 * prefix "instrument:complex:chip:cycle" and kUnknownRdbRunIdType.
 *
 * Note that some formats include the cycle_id and some do not.
 *
 * @param run_id Prefix input.
 * @return The formatted prefix as
 *         {date}:{instrument}:Q{queue_num}:R{run_num}:{cycle_id, optional per format}
 *         and the corresponding RdbRunIdType.
 */
std::tuple<std::string, RdbRunIdType> BuildNamePrefix(const std::string& run_id) {
  // HTP 1.1 format has 5 sections: date_instr_cd_group_num
  const auto parts = string::Split(run_id, "_");
  const auto num_parts = parts.size();
  if (kNumPartsHtp1OrGamma2 <= num_parts) {
    // this is the Gamma2 format, which is in the same order and format as the desired read output, but wihtout the
    // cycle
    if (parts[2].starts_with("Q") && parts[3].starts_with("R")) {
      return {fmt::format("{}:{}:{}:{}", parts[0], parts[1], parts[2], parts[3]), RdbRunIdType::kGamma2RdbRunIdType};
    }
    // this is the HTP format, date, group, run, instrument, chip, cycle
    //    220421_ENG-SYS-HTP_01_cloverleaf_WWX06R05C05_cycle01
    //    0      1           2  3          4            5
    //    date   group       │  instrument chip         cycle
    //                       run
    // The cycle field (parts[5]) is not included in the prefix; batch_index from the
    // RDB block header is used as the cycle_id by ParseRunIdAndCycleId, matching NSA behavior.
    return {fmt::format("{}:{}:{}:R{}", parts[0], parts[3], kUnspecifiedQueueNumber, parts[2]),
            RdbRunIdType::kHtp1RdbRunIdType};
  }
  if (parts.size() == kNumPartsDp2OrGamma1) {
    if (std::ranges::all_of(parts[4], [](const auto x) { return std::isdigit(x); })) {
      // If the 5th field is an integer, we have NS3_Gamma1 format.
      //    20230526_2-00063_13_123456A_1
      //    0        1       2  3       4
      //    date     |       |  tube    cycle
      //             |       run
      //             instrument
      // The bright-cycle (parts[4]) is not included in the prefix; batch_index from the
      // RDB block header is used as the cycle_id by ParseRunIdAndCycleId, matching NSA behavior.
      return {fmt::format("{}:{}:{}:R{}", parts[0], parts[1], kUnspecifiedQueueNumber, parts[2]),
              RdbRunIdType::kGamma1RdbRunIdType};
    }
    // If the 5th field is not an integer, we have DP2 format.
    //    220421_ENG-SYS-HTP_01_cloverleaf_WWX06R05C05
    //    0      1           2  3          4
    //    date   group       │  instrument chip
    //                       run
    return {fmt::format("{}:{}:{}:R{}", parts[0], parts[3], kUnspecifiedQueueNumber, parts[2]),
            RdbRunIdType::kDataPath2RdbRunIdType};
  }
  // using placeholder that adheres to the documented output format {date}:{instrument}:Q*:R*
  // cycle_id will be appended by LoadBlock() for kUnknownRdbRunIdType
  return {"date:instrument:Q*:R*", RdbRunIdType::kUnknownRdbRunIdType};
}

// NSA numeric base-64 alphabet – matches formatInteger<6, char*, uint32, 64, UNIFORM_UID_SIZE>.
// Order: 0-9, A-Z, a-z, +, _
static constexpr char kBase64Chars[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+_";

/**
 * Encode a u32 value to a base-64 string.
 *
 *   - Alphabet: 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+_
 *   - Pure numeric base-64 conversion (most-significant digit first)
 *   - Always produces exactly kMaxBase64UIDLength (6) characters, zero-padded
 *
 * 64^6 = 68,719,476,736 > 2^32, so every u32 value maps to a unique 6-character string.
 *
 * @param value The u32 cell index to encode
 * @return 6-character base-64 UID string
 */
std::string EncodeU32ToBase64(const u32 value) {
  std::string result(kMaxBase64UIDLength, '\0');
  u32 v = value;

  // Fill from right-to-left (least significant to most significant)
  // v & 0x3F is equivalent to v % 64
  // v >>= 6   is equivalent to v / 64
  for (s32 i = 5; i >= 0; --i) {
    result[i] = kBase64Chars[v & 0x3Fu];
    v >>= 6;
  }
  return result;
}

RdbSequenceReader::RdbSequenceReader(const fs::path& sequence_file_path)
    : _file_path(sequence_file_path),
      _read_seq_stream(std::make_unique<MemoryFileStream<>>(sequence_file_path.string())),
      _read_qual_stream(std::make_unique<MemoryFileStream<>>(sequence_file_path.string())),
      _read_len_stream(std::make_unique<MemoryFileStream<>>(sequence_file_path.string())),
      _read_cell_index_stream(std::make_unique<MemoryFileStream<>>(sequence_file_path.string())) {
  _file.exceptions(std::ifstream::failbit | std::ifstream::badbit | std::ifstream::eofbit);
  _file.open(sequence_file_path, std::ios::binary);

  _file.seekg(0, std::ios::end);
  _file_size_bytes = _file.tellg();
  if (_file_size_bytes <= 0) {
    // We are logging here in case there is a crash prior to the error from another thread
    const auto error_message = fmt::format("Error: '{}' is empty or malformed.", sequence_file_path.string());
    Logging::Error(error_message);
    throw error::Error(error_message);
  }

  _file.seekg(0, std::ios::beg);

  _read_seq_stream->SetFileSize(_file_size_bytes);
  _read_qual_stream->SetFileSize(_file_size_bytes);
  _read_len_stream->SetFileSize(_file_size_bytes);
  _read_cell_index_stream->SetFileSize(_file_size_bytes);
  _cur_block_start = 0;

  // load first block
  _cur_block_size = 0;
  _num_reads_in_block = 0;
  _cur_read_index = 0;
  _cur_batch_index = 0;
  LoadBlock();
}

// Validate that total_dataset_bytes and its rounded-up value fit within the remaining file.
static void ValidateBlockSize(const u64 total_dataset_bytes, const u64 remaining_bytes, const fs::path& file_path,
                              const std::streamoff block_start) {
  if (total_dataset_bytes > remaining_bytes) {
    throw error::Error(
        "Corrupted RDB block at offset {} in '{}': total_dataset_bytes ({}) exceeds remaining file size ({})",
        block_start, file_path.string(), total_dataset_bytes, remaining_bytes);
  }
  const auto max_streamoff_u64 = static_cast<u64>(std::numeric_limits<std::streamoff>::max());
  const auto rounded = RoundUp(total_dataset_bytes, static_cast<u64>(kChunkSize));
  if (rounded > remaining_bytes || rounded > max_streamoff_u64 - kChunkSize) {
    throw error::Error(
        "Corrupted RDB block at offset {} in '{}': rounded total_dataset_bytes ({}) exceeds remaining file size ({})",
        block_start, file_path.string(), rounded, remaining_bytes);
  }
}

/**
 * @brief Parse the run ID and cycle ID from the RDB block header.
 *
 * Extracts the run ID string from the block header, determines the RDB run ID format
 * (HTP 1.0, Datapath 2.0, NS3 Gamma1, NS3 Gamma2, or unknown) via @ref BuildNamePrefix,
 * and populates @c _read_name_prefix and @c _cycle_id accordingly:
 *
 * - **All types**: cycle ID is taken from the current batch index (@c _cur_batch_index)
 *   and appended to @c _read_name_prefix.  This matches NSA, which always uses the RDB
 *   block header's @c fullBrightCycleIndex as the cycle disambiguator.
 *
 * Also logs a one-time INFO message identifying the detected RDB type and app version.
 *
 * If the run ID string does not match any known format, @ref BuildNamePrefix falls back to
 * @c kUnknownRdbRunIdType with the placeholder prefix @c "date:instrument:Q*:R*"; this method
 * does not throw in that case.
 *
 * @param header A @c kChunkSize-byte view of the raw RDB block header.
 * @return A human-readable string identifying the RDB run ID type (e.g. @c "gamma2",
 *         @c "htp1", @c "gamma1", @c "data-path2", or @c "unknown"), suitable for logging.
 */
std::string RdbSequenceReader::ParseRunIdAndCycleId(const std::string_view header) {
  const auto app_version = ExtractHeaderString(header, kAppVersionOffset, kAppVersionLength);

  // Parse the run_id to identify the run type
  const auto* const run_id_start = header.data() + kRunIdOffset;
  const std::string run_id(run_id_start, std::find(run_id_start, run_id_start + kRunIdLength, '\0'));
  const auto [prefix, run_id_type] = BuildNamePrefix(run_id);

  // Format the run_id type name for logging (remove the -run-id-type suffix)
  auto run_id_formatted = enum_util::FormatEnumName(run_id_type);
  if (run_id_formatted.ends_with(kRunIdTypeSuffix)) {
    run_id_formatted.erase(run_id_formatted.size() - kRunIdTypeSuffix.size());
  }
  XOOS_LOG_INFO_ONCE("Detected RDB input of type '{}', app version '{}'", run_id_formatted, app_version);

  _read_name_prefix = prefix;

  // All RDB types use batch_index from the block header as the cycle disambiguator,
  // matching NSA behavior.  The legacy DataPath2 companion-file approach
  // (main_000001.rdb → expstate::cycle_id) produced a single value per directory,
  // causing duplicate read names when multiple RDB files share a directory.
  _cycle_id = _cur_batch_index;
  _read_name_prefix += fmt::format(":{}", _cycle_id.value());
  return run_id_formatted;
}

void RdbSequenceReader::LoadBlock() {
  char header[kChunkSize];

  try {
    _file.seekg(_cur_block_start, std::ios::beg);
    _file.read(header, kChunkSize);
  } catch (const std::ios_base::failure& e) {
    const auto error_message =
        fmt::format("Failed to read header from file '{}'.  Error: '{}'", _file_path.string(), e.what());
    // logging the error here in case there is a crash before the error message is displayed due to another thread's
    // failure
    Logging::Error(error_message);
    throw error::Error(error_message);
  }

  if (GetField<u16>(header, kBlockTypeOffset) != kNormalBlock) {
    throw error::Error("Not a normal block at {}", _cur_block_start);
  }

  const auto max_streamoff_u64 = static_cast<u64>(std::numeric_limits<std::streamoff>::max());
  if (_file_size_bytes < 0 || _cur_block_start < 0) {
    throw error::Error("Corrupted RDB reader state for '{}': negative file size ({}) or block offset ({})",
                       _file_path.string(), _file_size_bytes, _cur_block_start);
  }
  const auto file_size_u64 = static_cast<u64>(_file_size_bytes);
  const auto cur_block_start_u64 = static_cast<u64>(_cur_block_start);
  if (cur_block_start_u64 > file_size_u64 || file_size_u64 - cur_block_start_u64 < kChunkSize) {
    throw error::Error("Corrupted RDB block at offset {}: header exceeds remaining file size", _cur_block_start);
  }

  _num_reads_in_block = 0;
  _cur_read_index = 0;
  _cur_batch_index = GetField<u32>(header, kBatchIndexOffset);
  const auto header_view = std::string_view{header, kChunkSize};
  const auto run_id_formatted = ParseRunIdAndCycleId(header_view);

  const auto total_dataset_bytes = GetField<u64>(header, kTotalDatasetBytesOffset);
  if (total_dataset_bytes == 0) {
    _cur_block_size = kChunkSize;
    // no data in this block
    return;
  }

  const auto dataset_section_start_u64 = cur_block_start_u64 + kChunkSize;
  const auto remaining_dataset_bytes = file_size_u64 - dataset_section_start_u64;
  ValidateBlockSize(total_dataset_bytes, remaining_dataset_bytes, _file_path, _cur_block_start);
  const auto rounded_total_dataset_bytes = RoundUp(total_dataset_bytes, static_cast<u64>(kChunkSize));
  _cur_block_size = kChunkSize + static_cast<std::streamoff>(rounded_total_dataset_bytes);

  const auto num_datasets = GetField<u16>(header, kNumDatasetsOffset);
  if (num_datasets > kMaxNumDatasets) {
    throw error::Error("Corrupted RDB block at offset {}: num_datasets ({}) exceeds maximum ({})", _cur_block_start,
                       num_datasets, kMaxNumDatasets);
  }
  auto dataset_offset_u64 = dataset_section_start_u64;
  auto dataset_info = header + kBlockInfoSize;

  // Track which required datasets have been found by index into kRequiredDatasetNames.
  std::array<bool, kNumRequiredDatasets> dataset_found{};

  for (u16 i = 0; i < num_datasets; ++i, dataset_info += kDatasetInfoSize) {
    const auto dataset_name = dataset_info;
    if (dataset_offset_u64 > file_size_u64 || dataset_offset_u64 > max_streamoff_u64) {
      throw error::Error("Corrupted dataset offset in '{}': dataset '{}' starts at {} beyond file size {}",
                         _file_path.string(), dataset_name, dataset_offset_u64, file_size_u64);
    }
    const auto dataset_size_u64 = GetField<u64>(dataset_info, kDatasetSizeOffset);
    const auto remaining_bytes = file_size_u64 - dataset_offset_u64;
    if (dataset_size_u64 > remaining_bytes || dataset_size_u64 > max_streamoff_u64) {
      throw error::Error("Corrupted dataset size in '{}': dataset '{}' claims {} bytes but only {} remain",
                         _file_path.string(), dataset_name, dataset_size_u64, remaining_bytes);
    }
    const auto dataset_offset = static_cast<std::streamoff>(dataset_offset_u64);
    if (strcmp(dataset_name, "molecular_traces_count") == 0) {
      _file.seekg(dataset_offset, std::ios::beg);
      _file.read(reinterpret_cast<char*>(&_num_reads_in_block), sizeof(u32));
      dataset_found[kDatasetIdxTracesCount] = true;
    } else if (strcmp(dataset_name, "molecular_traces") == 0) {
      _read_seq_stream->Reset(dataset_offset);
      dataset_found[kDatasetIdxTraces] = true;
    } else if (strcmp(dataset_name, "molecular_traces_qscores") == 0) {
      _read_qual_stream->Reset(dataset_offset);
      dataset_found[kDatasetIdxQscores] = true;
    } else if (strcmp(dataset_name, "molecular_traces_read_length") == 0) {
      _read_len_stream->Reset(dataset_offset);
      dataset_found[kDatasetIdxReadLength] = true;
    } else if (strcmp(dataset_name, "mt_cell_index") == 0) {
      _read_cell_index_stream->Reset(dataset_offset);
      dataset_found[kDatasetIdxCellIndex] = true;
    }
    const auto rounded_dataset_size = RoundUp(dataset_size_u64, 8u);
    if (rounded_dataset_size > file_size_u64 - dataset_offset_u64) {
      throw error::Error("Corrupted dataset size in '{}': dataset '{}' rounded size {} exceeds {} remaining bytes",
                         _file_path.string(), dataset_name, rounded_dataset_size, file_size_u64 - dataset_offset_u64);
    }
    dataset_offset_u64 += rounded_dataset_size;
  }

  ValidateRequiredDatasets(dataset_found, _file_path, _cur_block_start, run_id_formatted, header_view);
}

BatchStatistics RdbSequenceReader::ReadBatchIntoArena(FixedReadRecordBatch& batch) {
  batch.start_time = std::chrono::high_resolution_clock::now();
  batch.num_bytes = 0ul;
  batch.num_bases = 0ul;
  for (batch.num_records = 0ul; batch.num_records < batch.Capacity() && HasMoreReads(); ++batch.num_records) {
    auto& record = (*batch.records)[batch.num_records];
    batch.num_bytes += ReadSingleFixed(record);
    batch.num_bases += record.SeqLen();
  }
  batch.end_time = std::chrono::high_resolution_clock::now();
  return {batch.num_records, batch.num_bases};
}

u32 RdbSequenceReader::ReadSingleFixed(FixedReadRecord& rec) {
  if (_cur_read_index >= _num_reads_in_block) {
    _cur_block_start += _cur_block_size;
    LoadBlock();
  }

  rec.Clear();

  const auto read_length = _read_len_stream->Read<u32>();
  // Guard: a zero-length read is never valid in practice and is a strong indicator of stream
  // desynchronization or file corruption.  Fail fast instead of silently producing empty records
  // that waste I/O cycles (which on networked storage can manifest as long CPU-near-0 stalls).
  if (read_length == 0) {
    auto error_string = fmt::format(
        "Read length of 0 encountered at read index {} in block at offset {} in file '{}'. "
        "This likely indicates file corruption or an incomplete RDB file.",
        _cur_read_index, _cur_block_start, _file_path.string());
    Logging::Error(error_string);
    throw error::Error(error_string);
  }
  auto read_length_packed = CeilDiv(read_length, 4u);
  u32 bytes_read = 0;

  if (read_length <= kMaxReadLength) {
    const auto read_cell_index = _read_cell_index_stream->Read<u32>();
    // We take the u32 and encode it to Base64. The buffer will always be there on the last 2 so we always exclude it by
    // truncating to 6 characters, providing a fixed length string for the read name.
    std::string b64_uid = EncodeU32ToBase64(read_cell_index);

    const auto [name_out, name_size] =
        fmt::format_to_n(rec.Name(), kBufferSize - 1, "{}:{}", _read_name_prefix, b64_uid);

    rec.comment_offset = name_out - rec.Name();

    rec.seq_offset = rec.comment_offset;
    _read_seq_stream->Read(rec.TwoBitsSeq(), read_length_packed);
    Reverse2BitOrder(rec.TwoBitsSeq(), read_length_packed);

    // fill padding regions with 0s
    std::memset(rec.TwoBitsSeq() - kTwoBitPadding, 0, kTwoBitPadding);
    // FIX: Was `+ read_length_packed + 1` which is off-by-one.  The +1 skips a byte and at
    // kMaxReadLength writes 1 byte past the end of the two_bit_seq[] array.
    // InitTwoBit() in read-record.h correctly uses `+ num_bytes` without +1.
    std::memset(rec.TwoBitsSeq() + read_length_packed, 0, kTwoBitPadding);

    DecodeDnaBases(rec.TwoBitsSeq(), rec.Seq(), read_length);

    rec.qual_offset = rec.seq_offset + read_length;
    _read_qual_stream->Read(_qual_buf, read_length_packed);
    DecodeQualScores(_qual_buf, rec.Qual(), read_length);

    rec.end_offset = rec.qual_offset + read_length;

    bytes_read = read_length_packed * 2;
    rec.SetStatus(FixedReadRecord::Status::kRead);
  } else {
    // FIX: Even though we discard this read, we must advance all per-read streams past its
    // data so that the next call to ReadSingleFixed() reads from the correct offsets.
    // Without these skips, every subsequent read in the block parses garbled data because
    // the cell-index, sequence, and quality streams are still pointing at *this* read's bytes.
    _read_cell_index_stream->Skip(sizeof(u32));
    _read_seq_stream->Skip(read_length_packed);
    _read_qual_stream->Skip(read_length_packed);
    rec.SetStatus(FixedReadRecord::Status::kTooLongFail);
  }
  ++_cur_read_index;
  return bytes_read;
}

bool RdbSequenceReader::HasMoreReads() const {
  return (_cur_read_index < _num_reads_in_block) || (_cur_block_start + _cur_block_size < _file_size_bytes);
}

}  // namespace xoos::demux
