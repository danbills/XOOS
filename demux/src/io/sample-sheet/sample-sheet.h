#pragma once

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "sequence/matcher/match-info.h"

namespace xoos::demux {

namespace fs = std::filesystem;

constexpr auto kSampleNameColumnHeader = "sample_name";
constexpr auto kSampleSidColumnHeader = "sample_sid";

/// @brief (SID sequence, sample name) pairs in CSV/TSV row order.
using SidSeqNamePairs = std::vector<std::pair<std::string, std::string>>;

/**
 * @brief Read a sample sheet (CSV or TSV) and return (sequence, sample_name) pairs in row order.
 *
 * Parses the file expecting columns named "sample_sid" and "sample_name". Validates that
 * neither SID sequences nor sample names are duplicated, and that no sample uses the
 * reserved name for unassigned reads.
 *
 * @param path Path to the sample sheet file.
 * @return Pairs of (SID sequence, sample name) preserving CSV/TSV row order.
 * @throws error::Error On duplicate SIDs, duplicate sample names, or reserved name usage.
 */
SidSeqNamePairs ReadSampleSheet(const fs::path& path);

/**
 * @brief Loads a sample sheet and converts it into a BarcodePool with contiguous 0-based IDs.
 *
 * Each row in the sample sheet (sequence, sample name) becomes a Barcode whose ID is assigned
 * sequentially in CSV row order. The resulting pool is passed to LoadLutBundle so that only the
 * user-specified SIDs are built into the LUT.
 *
 * @param sample_sheet Path to the sample sheet CSV/TSV file.
 * @return A BarcodePool with one entry per sample sheet row, IDs 0..N-1.
 */
BarcodePool LoadSampleSheet(const fs::path& sample_sheet);

}  // namespace xoos::demux
