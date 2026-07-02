#include "io/write-bigwig.h"

#include <fstream>
#include <sstream>
#include <vector>

#include <xoos/log/logging.h>

#include "observations.h"

extern "C" {
#include <libbigwig/bigWig.h>
}

namespace xoos::cnc {

constexpr int kBigWigInitBufferSize = 1 << 17;  // 128 KiB

std::pair<std::map<std::string, u32>, std::vector<std::string>> GetChromLengthsFromFaiFile(
    const std::string& fai_filename) {
  std::map<std::string, u32> chrom_lengths;
  std::vector<std::string> chrom_order;  // Preserve .fai order
  std::ifstream fai_file(fai_filename);

  if (!fai_file.is_open()) {
    throw std::runtime_error("Could not open .fai file: " + fai_filename);
  }

  std::string line;
  while (std::getline(fai_file, line)) {
    std::istringstream iss(line);
    std::string chrom_name;
    u32 length;

    if (iss >> chrom_name >> length) {
      chrom_lengths[chrom_name] = length;
      chrom_order.push_back(chrom_name);  // Preserve original order
    }
  }

  return {chrom_lengths, chrom_order};
}

bool ValidateIntervalsSorted(const std::vector<const char*>& chroms,
                             const std::vector<u32>& starts,
                             const std::vector<std::string>& chrom_order) {
  if (chroms.empty()) {
    return true;
  }

  // Create chromosome order map using the ACTUAL .fai order
  std::map<std::string, s32> chrom_order_map;
  for (size_t i = 0; i < chrom_order.size(); ++i) {
    chrom_order_map[chrom_order[i]] = static_cast<s32>(i);
  }

  for (size_t i = 1; i < chroms.size(); ++i) {
    std::string prev_chrom(chroms[i - 1]);
    std::string curr_chrom(chroms[i]);

    // Get chromosome order indices
    auto prev_it = chrom_order_map.find(prev_chrom);
    auto curr_it = chrom_order_map.find(curr_chrom);
    if (prev_it == chrom_order_map.end() || curr_it == chrom_order_map.end()) {
      Logging::Error("Unknown chromosome found: {} or {}", prev_chrom, curr_chrom);
      return false;
    }

    s32 prev_order = prev_it->second;
    s32 curr_order = curr_it->second;

    if (prev_order > curr_order) {
      Logging::Error("Chromosomes out of order at interval {}: {} (order {}) comes after {} (order {})",
                     i,
                     curr_chrom,
                     curr_order,
                     prev_chrom,
                     prev_order);
      return false;
    } else if (prev_order == curr_order) {
      // Same chromosome - check position order
      if (starts[i] < starts[i - 1]) {
        Logging::Error(
            "Positions out of order at interval {} on {}: {} comes after {}", i, curr_chrom, starts[i], starts[i - 1]);
        return false;
      }
    }
  }

  return true;
}

void WriteObservationsToBigWig(const Observations& obs, const fs::path& filename, const fs::path& fai_filename) {
  // Read chromosome information from .fai file
  std::map<std::string, u32> chrom_lengths;
  std::vector<std::string> chrom_order;
  try {
    auto [cl, co] = GetChromLengthsFromFaiFile(fai_filename.string());
    chrom_lengths = std::move(cl);
    chrom_order = std::move(co);
    Logging::Info("Read {} chromosomes from {}", chrom_lengths.size(), fai_filename.string());
  } catch (const std::exception& e) {
    throw std::runtime_error("Error reading .fai file: " + std::string(e.what()));
  }

  // Open BigWig file for writing
  bigWigFile_t* fp = bwOpen(filename.c_str(), nullptr, "w");
  if (fp == nullptr) {
    throw std::runtime_error("Failed to open BigWig file: " + filename.string());
  }

  try {
    // Create stable C-style arrays from .fai data using the CORRECT ORDER
    std::vector<std::string> chrom_name_strings;
    chrom_name_strings.reserve(chrom_order.size());  // Use chrom_order size

    // Build chrom_name_strings using chrom_order instead of iterating over the map
    for (const std::string& chrom_name : chrom_order) {
      chrom_name_strings.push_back(chrom_name);
    }

    // Create C arrays that won't be invalidated
    std::vector<const char*> chrom_names_c;
    std::vector<u32> chrom_lens_c;
    for (const auto& chrom_name : chrom_name_strings) {
      chrom_names_c.push_back(chrom_name.c_str());
      chrom_lens_c.push_back(chrom_lengths[chrom_name]);
    }
    // Create chromosome list and assign to file structure
    fp->cl = bwCreateChromList(chrom_names_c.data(), chrom_lens_c.data(), static_cast<s64>(chrom_names_c.size()));
    if (fp->cl == nullptr) {
      throw std::runtime_error("Failed to create chromosome list");
    }

    // Create and write header
    // Zoom level of 7 is chosen based on https://academic.oup.com/bib/article/14/2/178/208453
    if (bwCreateHdr(fp, 7) != 0) {
      throw std::runtime_error("Failed to create BigWig header");
    }

    if (bwWriteHdr(fp) != 0) {
      throw std::runtime_error("Failed to write BigWig header");
    }

    // Prepare data arrays from observations
    // FIRST: Create all the strings (reserve space to avoid reallocation)
    std::vector<std::string> chrom_strings;
    chrom_strings.reserve(obs.contigs.size());  // Prevent reallocation
    for (const auto& contig : obs.contigs) {
      chrom_strings.push_back(contig);
    }

    // SECOND: Create pointer arrays after all strings are stable
    std::vector<const char*> chroms;
    std::vector<u32> starts_vec;
    std::vector<u32> ends_vec;
    std::vector<float> values_vec;

    for (size_t i = 0; i < obs.contigs.size(); ++i) {
      // Check if chromosome exists in the reference
      if (chrom_lengths.find(chrom_strings[i]) == chrom_lengths.end()) {
        throw std::runtime_error("Chromosome '" + chrom_strings[i] + "' not found in .fai file");
      }
      chroms.push_back(chrom_strings[i].c_str());
      starts_vec.push_back(static_cast<u32>(obs.starts[i]));
      ends_vec.push_back(static_cast<u32>(obs.ends[i]));
      values_vec.push_back(static_cast<float>(obs.obvs[i]));
    }

    // Validate sorting before writing
    Logging::Info("Validating interval sorting for bigWig writing...");
    if (!ValidateIntervalsSorted(chroms, starts_vec, chrom_order)) {
      throw std::runtime_error("Intervals are not properly sorted for BigWig format");
    }

    // Write intervals to BigWig
    if (bwAddIntervals(fp, chroms.data(), starts_vec.data(), ends_vec.data(), values_vec.data(), obs.contigs.size()) !=
        0) {
      throw std::runtime_error("Failed to add intervals to BigWig");
    }

    // Close the file
    bwClose(fp);
    Logging::Info("BigWig file '{}' created successfully with {} intervals", filename.string(), obs.contigs.size());
  } catch (...) {
    // Ensure file is closed on any exception
    bwClose(fp);
    throw;
  }
}

void WriteEmptyBigWig(const fs::path& filename, const fs::path& fai_filename) {
  if (bwInit(kBigWigInitBufferSize) != 0) {
    throw std::runtime_error("Failed to initialize libbigwig");
  }

  auto* fp = bwOpen(filename.c_str(), nullptr, "w");
  if (fp == nullptr) {
    bwCleanup();
    throw std::runtime_error("Failed to open BigWig file for writing: " + filename.string());
  }

  try {
    auto [chrom_lengths, chrom_order] = GetChromLengthsFromFaiFile(fai_filename.string());

    std::vector<const char*> chrom_names_c;
    std::vector<u32> chrom_lens_c;
    for (const auto& name : chrom_order) {
      chrom_names_c.push_back(name.c_str());
      chrom_lens_c.push_back(chrom_lengths[name]);
    }

    fp->cl = bwCreateChromList(chrom_names_c.data(), chrom_lens_c.data(), static_cast<s64>(chrom_names_c.size()));
    if (fp->cl == nullptr) {
      throw std::runtime_error("Failed to create chromosome list for empty BigWig");
    }

    if (bwCreateHdr(fp, 0) != 0) {
      throw std::runtime_error("Failed to create BigWig header");
    }
    if (bwWriteHdr(fp) != 0) {
      throw std::runtime_error("Failed to write BigWig header");
    }

    bwClose(fp);
    bwCleanup();
  } catch (...) {
    bwClose(fp);
    bwCleanup();
    throw;
  }
}

}  // namespace xoos::cnc
