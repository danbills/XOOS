#pragma once

#include <filesystem>

namespace xoos {

/// RAII guard that removes a file on destruction.
/// If the file does not exist at construction time, it is created as an empty file.
class TmpFile {
 public:
  explicit TmpFile(std::filesystem::path path);

  ~TmpFile();

  TmpFile(const TmpFile&) = delete;
  TmpFile& operator=(const TmpFile&) = delete;
  TmpFile(TmpFile&&) = delete;
  TmpFile& operator=(TmpFile&&) = delete;

  [[nodiscard]] const std::filesystem::path& Path() const;

  /// Create an empty file at @p path, including any missing parent directories.
  static void CreateTmpFile(const std::filesystem::path& path);

 private:
  void Cleanup() const;

  std::filesystem::path _path;
};

}  // namespace xoos
