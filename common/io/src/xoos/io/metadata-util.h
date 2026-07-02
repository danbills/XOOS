#pragma once
#include <optional>
#include <string>

namespace xoos::io {

const std::string kTsvCommentLinePrefix = "#";
const std::string kVcfCommentLinePrefix = "##";
const std::string kRocheCommandLinePrefix = "##RocheCommandLine=<";
const std::string kKeyId = "ID";
const std::string kKeyVersion = "Version";
const std::string kKeyCommandLine = "CommandLine";

/**
 * @brief Holds program name, version, and full command line for metadata output.
 */
struct CommandLineInfo {
  std::string name;
  std::string version;
  std::string command_line;
};

/**
 * @brief Format a CommandLineInfo as a RocheCommandLine metadata line.
 *
 * Produces: ##RocheCommandLine=<ID=name,Version="version",CommandLine="command_line">
 *
 * Used for both VCF and TSV output.
 */
std::string GetCommandInfo(const CommandLineInfo& info);

/**
 * @brief Parse a RocheCommandLine metadata line back into a CommandLineInfo struct.
 *
 * Expected format:
 *   ##RocheCommandLine=<ID=name,Version="version",CommandLine="command_line">
 *
 * Key-value pairs inside the angle brackets are parsed in any order.
 * Handles escaped quotes and backslashes in quoted values (inverse of EscapeMetadataValue).
 * Returns std::nullopt silently if the line is not a RocheCommandLine line, or with a
 * warning log if the line starts with ##RocheCommandLine= but is malformed.
 */
std::optional<CommandLineInfo> ParseCommandInfo(const std::string& line);

/**
 * @brief Write a RocheCommandLine metadata comment to an output stream.
 *
 * Writes a single ##RocheCommandLine=<...> comment line directly to the stream
 * to avoid CSV quoting of special characters (commas, double-quotes).
 */
template <typename OutputStream>
void WriteTsvMetadata(OutputStream& stream, const CommandLineInfo& info) {
  stream << GetCommandInfo(info) << "\n";
}

}  // namespace xoos::io
