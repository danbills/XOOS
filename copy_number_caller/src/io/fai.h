#pragma once
#include <string>
#include <unordered_map>
#include <vector>

#include <xoos/types/fs.h>

namespace xoos::cnc {

using SeqLenMap = std::unordered_map<std::string, size_t>;
SeqLenMap ParseFai(const fs::path& fname);
std::unordered_map<std::string, size_t> GetContigOrder(const fs::path& fname);
std::vector<std::string> GetContigsInOrder(const fs::path& fname);
}  // namespace xoos::cnc
