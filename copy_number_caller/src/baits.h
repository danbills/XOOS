#pragma once
#include <armadillo>
#include <string>
#include <utility>
#include <vector>

#include <xoos/io/fasta-reader.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>

#include "bigwig_cpp.h"
#include "io/fai.h"

namespace xoos::cnc {
const f64 kCoverageMinOnTargetMapp = 0.6;
const f64 kCoverageMinOffTargetMapp = 0.1;
const f64 kCoverageGCQuantile = 0.001;

/**
 * @class BaitRecords
 * @brief container for records in a bait
 *
 * Input file should have the following format:
 * line 1: header
 * line 2-n: <region>\t<gc bias>\t<mappability>\t<rep timing>\t<Gene>\t<on_target>\n
 */
class BaitRecords {
 public:
  BaitRecords() = default;
  BaitRecords(const BaitRecords&) = default;
  BaitRecords(BaitRecords&&) noexcept = default;
  BaitRecords& operator=(const BaitRecords&) = default;
  BaitRecords& operator=(BaitRecords&&) noexcept = default;

  explicit BaitRecords(std::istream& ifs);

  const std::vector<std::string>& GetRegions() const {
    return _region;
  }

  const arma::vec& GetMappability() const {
    return _mappability;
  }

  const arma::vec& GetGCBias() const {
    return _gc_bias;
  }

  const std::vector<bool>& GetOnTargetStatus() const {
    return _on_target;
  }

  void SetRegions(std::vector<std::string>&& regions) {
    _region = std::move(regions);
  }

  void SetMappability(const arma::vec& mappability) {
    _mappability = mappability;
  }

  void SetMappability(arma::vec&& mappability) {
    _mappability = std::move(mappability);
  }

  void SetGCBias(const arma::vec& gc_bias) {
    _gc_bias = gc_bias;
  }

  void SetGCBias(arma::vec&& gc_bias) {
    _gc_bias = std::move(gc_bias);
  }

  void SetOnTargetStatus(const std::vector<bool>& on_target) {
    _on_target = on_target;
  }

  void SetOnTargetStatus(std::vector<bool>&& on_target) {
    _on_target = std::move(on_target);
  }

  void SetReferenceFile(const std::string& reference_name) {
    _reference_name = reference_name;
  }

  void SetAllOnTarget();

  void SetSeqLengths(SeqLenMap& seq_lengths) {
    _seq_lengths = seq_lengths;
  }

  void SetSeqLengthsFromHeader(const std::string& header);

  void SetReferenceSequenceFromHeader(const std::string& header);

  const std::string& GetReferenceFile() const {
    return _reference_name;
  }

  std::vector<std::string> GetOnOrOffTargetRegions(bool);

  const SeqLenMap& GetSeqLengths() const {
    return _seq_lengths;
  }

  void SortByRegion(const fs::path& fai_file);
  void AnnotateGC(io::FastaReader& fa);
  void AnnotateMappability(BigWig& bigwig);

  void Write(std::ofstream& ofs, const io::CommandLineInfo& command_line_info) const;

  // Remove baits with GC content outside the specified thresholds
  void FilterExtremeGC(f64 lower_threshold = 0.15, f64 upper_threshold = 0.75);

 private:
  void LoadBaitFile(std::istream& ifs);
  std::vector<std::string> _region{};
  arma::vec _gc_bias;
  arma::vec _mappability;
  std::vector<bool> _on_target{};
  SeqLenMap _seq_lengths{};
  std::string _reference_name;
};

}  // namespace xoos::cnc
