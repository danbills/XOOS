#include "xoos/io/metadata-util.h"

#include <fmt/format.h>

#include <xoos/log/logging.h>

namespace xoos::io {

static std::string EscapeMetadataValue(const std::string& value) {
  std::string escaped;
  escaped.reserve(value.size());
  for (char c : value) {
    if (c == '\\' || c == '"') {
      escaped += '\\';
    }
    escaped += c;
  }
  return escaped;
}

std::string GetCommandInfo(const CommandLineInfo& info) {
  return fmt::format(R"(##RocheCommandLine=<ID={},Version="{}",CommandLine="{}">)",
                     info.name,
                     EscapeMetadataValue(info.version),
                     EscapeMetadataValue(info.command_line));
}

/**
 * @brief Unescape a value that was escaped by EscapeMetadataValue.
 *
 * Reverses backslash-escaping of '"' and '\' characters.
 */
static std::string UnescapeMetadataValue(const std::string& value) {
  std::string result;
  result.reserve(value.size());
  for (size_t i = 0; i < value.size(); ++i) {
    if (value[i] == '\\' && i + 1 < value.size()) {
      result += value[i + 1];
      ++i;
    } else {
      result += value[i];
    }
  }
  return result;
}

/**
 * @brief Strip surrounding double quotes from a value and unescape the interior.
 *
 * Returns the unescaped content between the quotes, or the original value if it is not quoted.
 */
static std::string StripQuotes(const std::string& value) {
  if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
    return UnescapeMetadataValue(value.substr(1, value.size() - 2));
  }
  return value;
}

/**
 * @brief Split the body of a RocheCommandLine tag by commas, respecting quoted values.
 *
 * Commas inside double-quoted strings (including escaped quotes) are not treated as delimiters.
 */
static std::vector<std::string> SplitKeyValuePairs(const std::string& body) {
  std::vector<std::string> tokens;
  std::string current;
  bool in_quotes = false;
  for (size_t i = 0; i < body.size(); ++i) {
    if (body[i] == '\\' && in_quotes && i + 1 < body.size()) {
      current += body[i];
      current += body[i + 1];
      ++i;
    } else if (body[i] == '"') {
      in_quotes = !in_quotes;
      current += body[i];
    } else if (body[i] == ',' && !in_quotes) {
      tokens.push_back(current);
      current.clear();
    } else {
      current += body[i];
    }
  }
  if (!current.empty()) {
    tokens.push_back(current);
  }
  return tokens;
}

std::optional<CommandLineInfo> ParseCommandInfo(const std::string& line) {
  if (line.rfind(kRocheCommandLinePrefix, 0) != 0) {
    return std::nullopt;
  }
  if (line.empty() || line.back() != '>') {
    Logging::Warn("Malformed RocheCommandLine line: missing closing '>'");
    return std::nullopt;
  }

  auto body = line.substr(kRocheCommandLinePrefix.size(), line.size() - kRocheCommandLinePrefix.size() - 1);
  auto tokens = SplitKeyValuePairs(body);

  CommandLineInfo info;
  bool has_id = false;
  bool has_version = false;
  bool has_command_line = false;

  for (const auto& token : tokens) {
    auto eq_pos = token.find('=');
    if (eq_pos == std::string::npos) {
      Logging::Warn("Malformed RocheCommandLine token (no '='): '{}'", token);
      return std::nullopt;
    }
    auto key = token.substr(0, eq_pos);
    auto value = token.substr(eq_pos + 1);

    if (key == kKeyId) {
      info.name = value;
      has_id = true;
    } else if (key == kKeyVersion) {
      info.version = StripQuotes(value);
      has_version = true;
    } else if (key == kKeyCommandLine) {
      info.command_line = StripQuotes(value);
      has_command_line = true;
    }
  }

  if (!has_id || !has_version || !has_command_line) {
    Logging::Warn("Malformed RocheCommandLine line: missing required field(s) (ID, Version, CommandLine)");
    return std::nullopt;
  }

  return info;
}

}  // namespace xoos::io
