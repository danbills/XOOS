#include "baits.h"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include <fmt/core.h>

#include <csv.hpp>

#include <xoos/error/error.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/util/container-functions.h>
#include <xoos/util/parse-int.h>

#include "io/column-names.h"
#include "utility/utility-functions.h"

namespace xoos::cnc {

/**
 *  Bait file format described in README.md
 */
BaitRecords::BaitRecords(std::istream& ifs) {
  LoadBaitFile(ifs);
}

/**
 * @brief loads a baits BED file.
 * Expected BED format: Contig\tStart\tEnd\tGCBias\tMappability\tOnTarget (0-based half-open).
 * `#`-prefixed metadata/comments and a `#`-prefixed column header are supported and skipped.
 * @param ifs input stream
 */
void BaitRecords::LoadBaitFile(std::istream& ifs) {
  std::vector<f64> gc_biases;
  std::vector<f64> mappabilities;
  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::string contig;
    std::string start_str;
    std::string end_str;
    std::string gc_str;
    std::string mapp_str;
    std::string on_target_str;
    std::istringstream ss(line);
    if (!std::getline(ss, contig, '\t') || !std::getline(ss, start_str, '\t') || !std::getline(ss, end_str, '\t') ||
        !std::getline(ss, gc_str, '\t') || !std::getline(ss, mapp_str, '\t') ||
        !std::getline(ss, on_target_str, '\t')) {
      throw error::Error(
          "Baits BED file must have 6 columns (Contig, Start, End, GCBias, Mappability, OnTarget): \"{}\"", line);
    }
    // reconstruct 1-based region string for internal use
    arma::uword start, end;
    try {
      start = std::stoul(start_str);
      end = std::stoul(end_str);
    } catch (const std::exception& e) {
      throw error::Error("Baits BED file parse error: {} (line: \"{}\")", e.what(), line);
    }
    std::string r = contig + ":" + std::to_string(start + 1) + "-" + std::to_string(end);
    _region.push_back(r);
    gc_biases.push_back(std::stod(gc_str));
    mappabilities.push_back(std::stod(mapp_str));
    _on_target.push_back(on_target_str == "TRUE");
  }
  _gc_bias = arma::vec(gc_biases);
  _mappability = arma::vec(mappabilities);
}

/**
 * @brief set the SeqLengths member, given a SAM header with @SQ lines (see SAM format or GATK interval format)
 * @param header
 */
void BaitRecords::SetSeqLengthsFromHeader(const std::string& header) {
  std::istringstream is(header);
  std::string line;
  while (std::getline(is, line)) {
    if (line[0] != '@') {
      throw std::runtime_error("header has a line not prepended with '@'!");
    } else {
      if (line.substr(0, 3) == "@SQ") {
        std::istringstream line_is(line);
        std::string field;
        std::string contig;
        size_t length = 0;
        while (line_is) {
          line_is >> field;
          if (field.substr(0, 3) == "SN:") {
            contig = field.substr(3, field.size() - 3);
          } else if (field.substr(0, 3) == "LN:") {
            length = util::ParseU64(field.substr(3, field.size() - 3));
          }
        }
        if (!contig.empty() && length != 0) {
          _seq_lengths[contig] = length;
        } else {
          throw error::Error("invalid contig and/or length in header");
        }
      }
    }
  }
}

void BaitRecords::SetReferenceSequenceFromHeader(const std::string& header) {
  size_t ur_pos = header.find("UR:");
  if (ur_pos == std::string::npos) {
    Logging::Info("No reference sequence UR tag found in header");
    _reference_name = "";
  } else {
    std::string tmp = header.substr(ur_pos);
    size_t whitespace_pos = tmp.find_first_of("\t\n");
    _reference_name = tmp.substr(3, whitespace_pos - 3);
    Logging::Info("found UR: {}", _reference_name);
  }
}

/**
 * @brief get choice of on-target OR off-target intervals
 * @param o if true, return on-targets, if false return off-targets
 * @return a BaitRecords object containing the desired target intetvals
 */
std::vector<std::string> BaitRecords::GetOnOrOffTargetRegions(bool o) {
  std::vector<std::string> ret;
  for (size_t i = 0; i < _region.size(); ++i) {
    if ((o && _on_target[i]) || (!o && !_on_target[i])) {
      ret.push_back(_region[i]);
    }
  }
  return ret;
}

/**
 * @brief set all targets to have on-target status
 */
void BaitRecords::SetAllOnTarget() {
  _on_target = std::vector<bool>(_region.size(), true);
}

/**
 * @brief sort the object by contig and position, contig order defined by a provided FAI file
 * @param fai_file fai file containing desired order of contigs
 */
void BaitRecords::SortByRegion(const fs::path& fai_file) {
  if (_region.empty()) {
    Logging::Error("cannot SortByRegion because there are no regions for this BaitRecords object!");
    throw std::runtime_error("empty region std::vector");
  }
  std::unordered_map<std::string, size_t> contig_order = GetContigOrder(fai_file);
  // this will contain the "correct" order by which to reorder everythign else
  std::vector<arma::uword> argsort(_region.size());
  std::iota(argsort.begin(), argsort.end(), 0);
  std::stable_sort(argsort.begin(), argsort.end(), [this, &contig_order](size_t i, size_t j) {
    const auto& [lcontig, lstart, lend] = ParseRegionString(_region[i]);
    const auto& [rcontig, rstart, rend] = ParseRegionString(_region[j]);
    return std::tie(contig_order[lcontig], lstart) < std::tie(contig_order[rcontig], rstart);
  });
  // sort all the member std::vectors
  _region = util::container::Permute(_region, argsort);
  if (!_on_target.empty()) {
    _on_target = util::container::Permute(_on_target, argsort);
  }
  arma::uvec argsort_arma(argsort);
  if (!_gc_bias.empty()) {
    _gc_bias = _gc_bias.elem(argsort_arma);
  }
  if (!_mappability.empty()) {
    _mappability.elem(argsort_arma);
  }
}

/**
 * @brief annotate intervals with their respective GC contents
 * @param fa a FastaReader object
 */
void BaitRecords::AnnotateGC(io::FastaReader& fa) {
  _gc_bias.clear();
  _gc_bias.resize(_region.size());
  for (size_t i = 0; i < _region.size(); ++i) {
    const auto& [contig, start, end] = ParseRegionString(_region[i]);
    std::string seq(fa.GetSequence(contig, static_cast<s32>(start), static_cast<s32>(end), true));
    size_t gc_count = std::ranges::count_if(seq.begin(), seq.end(), [](char c) { return c == 'G' || c == 'C'; });
    _gc_bias[i] = static_cast<f64>(gc_count) / static_cast<f64>(seq.size());
  }
}

void BaitRecords::FilterExtremeGC(f64 lower_threshold, f64 upper_threshold) {
  if (_gc_bias.empty()) {
    throw std::runtime_error("GC content not annotated. Call AnnotateGC() first.");
  }

  size_t original_size = _region.size();

  // Filter by creating new std::vectors with only kept elements
  std::vector<std::string> filtered_region;
  std::vector<bool> filtered_on_target;
  std::vector<f64> filtered_gc;
  std::vector<f64> filtered_mapp;

  for (size_t i = 0; i < _region.size(); ++i) {
    if (_gc_bias[i] >= lower_threshold && _gc_bias[i] <= upper_threshold) {
      filtered_region.push_back(_region[i]);
      filtered_on_target.push_back(_on_target[i]);
      filtered_gc.push_back(_gc_bias[i]);
      if (!_mappability.empty()) {
        filtered_mapp.push_back(_mappability[i]);
      }
    }
  }

  // Replace with filtered std::vectors
  _region = std::move(filtered_region);
  _on_target = std::move(filtered_on_target);
  _gc_bias = arma::vec(filtered_gc);
  if (!_mappability.empty()) {
    _mappability = arma::vec(filtered_mapp);
  }

  size_t removed = original_size - _region.size();
  Logging::Info("Removed {} baits with extreme GC content (kept GC range {:.2f}-{:.2f})",
                removed,
                lower_threshold,
                upper_threshold);
}

/**
 * @brief annotate intervals with their respective mappabilities
 * @param mapp_bigwig
 */
void BaitRecords::AnnotateMappability(BigWig& mapp_bigwig) {
  _mappability.clear();
  _mappability.resize(_region.size());
  for (size_t i = 0; i < _region.size(); ++i) {
    const auto& [contig, start, end] = ParseRegionString(_region[i]);
    _mappability[i] = mapp_bigwig.GetMean(contig, start, end);
  }
}

/**
 * @brief write BaitRecords to BED file with metadata and column header.
 * @param ofs output file stream
 * @param command_line_info command line metadata including program version and rendered command line
 */
void BaitRecords::Write(std::ofstream& ofs, const io::CommandLineInfo& command_line_info) const {
  auto writer = csv::make_tsv_writer(ofs);
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info);
  }
  vec<std::string> column_names{
      kColumnContigAsComment, kColumnStart, kColumnEnd, kColumnGCBias, kColumnMappability, kColumnOnTarget};
  writer << column_names;
  for (size_t i = 0; i < _region.size(); ++i) {
    const auto& [contig, start, end] = ParseRegionString(_region[i]);
    writer << std::make_tuple(contig, start, end, _gc_bias[i], _mappability[i], _on_target[i] ? "TRUE" : "FALSE");
  }
}

}  // namespace xoos::cnc
