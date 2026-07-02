#pragma once
#include <map>
#include <string>
#include <vector>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>

extern "C" {
#include <libbigwig/bigWig.h>
}

namespace xoos::cnc {

// Forward declaration
class Observations;

/**
 * @brief Reads chromosome lengths from a FASTA index (.fai) file
 *
 * This function parses a FASTA index file to extract chromosome names and their
 * corresponding lengths. The .fai file format contains tab-separated values where
 * the first column is the chromosome name and the second column is its length.
 *
 * @param fai_filename fs::path to the FASTA index (.fai) file to read
 * @return std::pair containing:
 *         - first: A map with chromosome names as keys and their lengths as values
 *         - second: A std::vector of chromosome names in the order they appear in the file
 * @throws std::runtime_error if the file cannot be opened or parsed
 */
std::pair<std::map<std::string, u32>, std::vector<std::string>> GetChromLengthsFromFaiFile(
    const std::string& fai_filename);

/**
 * @brief Validates that genomic intervals are sorted according to chromosome order and genomic positions.
 *
 * This function checks whether a collection of genomic intervals is properly sorted first by
 * chromosome order (as specified in the reference genome index) and then by genomic position
 * within each chromosome. The validation ensures that intervals can be processed sequentially
 * for formats that require sorted input (e.g., BigWig files).
 *
 * @param chroms Vector of chromosome names as C-strings for each interval
 * @param starts Vector of start positions (0-based) corresponding to each interval
 * @param chrom_order Vector defining the canonical chromosome ordering (typically from .fai file)
 *
 * @return true if intervals are properly sorted, false otherwise
 *
 * @details The function performs the following validation steps:
 * - Returns true immediately for empty input
 * - Creates a mapping from chromosome names to their canonical order indices
 * - Iterates through consecutive interval pairs checking:
 *   - Both chromosomes exist in the reference order
 *   - Current chromosome order >= previous chromosome order
 *   - For same chromosome: current position >= previous position
 * - Logs detailed error messages for any ordering violations
 *
 * @note This function assumes chroms and starts std::vectors have the same size
 * @warning Will return false and log errors if any chromosome is not found in chrom_order
 */
bool ValidateIntervalsSorted(const std::vector<const char*>& chroms,
                             const std::vector<u32>& starts,
                             const std::vector<std::string>& chrom_order);

/**
 * @brief Write the observations to a BigWig file
 * @param obs The observations object to write
 * @param filename path to output BigWig file
 * @param fai_filename path to FASTA index file containing chromosome lengths
 */
void WriteObservationsToBigWig(const Observations& obs, const fs::path& filename, const fs::path& fai_filename);

/**
 * @brief Write a valid but data-empty BigWig file (header + chromosome list, no intervals).
 *
 * Used by the low-coverage early-exit path to produce structurally valid BigWig
 * files that downstream consumers (IGV, pyBigWig) can open without error.
 *
 * @param filename path to output BigWig file
 * @param fai_filename path to FASTA index file containing chromosome lengths
 */
void WriteEmptyBigWig(const fs::path& filename, const fs::path& fai_filename);

}  // namespace xoos::cnc
