#pragma once

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <htslib/vcf.h>

#include <xoos/io/htslib-util/kstring.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "vcf-header.h"

namespace xoos::io {
constexpr auto kGT = "GT";

inline bool IsGtStandardFormatField(const std::string& field) {
  return field == kGT;
}

class VcfRecord;
using VcfRecordPtr = std::shared_ptr<VcfRecord>;
using BcfRecordPtr = std::shared_ptr<bcf1_t>;

class VcfRecord {
 public:
  explicit VcfRecord(const VcfHeaderPtr& hdr);
  static VcfRecordPtr CreateFromHeader(const VcfHeaderPtr& hdr, s32 unpack);
  static VcfRecordPtr CreateFromHeader(const VcfHeaderPtr& hdr);
  static VcfRecordPtr ReadFromFile(const VcfHeaderPtr& hdr, const HtsFileSharedPtr& input_vcf_name, s32 unpack);
  static VcfRecordPtr ReadFromFile(const VcfHeaderPtr& hdr, const HtsFileSharedPtr& input_vcf_name);
  static VcfRecordPtr ReadFromRegion(const VcfHeaderPtr& hdr,
                                     const HtsFileSharedPtr& input_vcf_fp,
                                     tbx_t* tbx_idx,
                                     hts_itr_t* hts_itr,
                                     Kstring& kstr,
                                     s32 unpack);

  std::string Chromosome() const;
  hts_pos_t Position() const;
  bool IsSnp() const;
  bool IsPass() const;
  std::string Id() const;
  std::string Allele(s32 which) const;
  s32 NumAlleles() const;

  bool HasInfoFieldNoCheck(const std::string& field) const;
  template <typename T>
  std::vector<T> GetInfoFieldNoCheck(const std::string& field) const;
  std::string GetInfoFieldStringNoCheck(const std::string& field) const;
  bool HasFormatFieldNoCheck(const std::string& field) const;
  template <typename T>
  std::vector<T> GetFormatFieldNoCheck(const std::string& field) const;
  std::string GetFormatFieldStringNoCheck(const std::string& field) const;

  void ConfirmFieldInHeader(const std::string& field, const std::string& field_type) const;
  void ConfirmInfoFieldInHeader(const std::string& field) const;
  void ConfirmFormatFieldInHeader(const std::string& field) const;

  template <typename T>
  std::vector<T> GetInfoField(const std::string& field) const;
  void SetInfoField(const std::string& field, const std::vector<f32>& values);
  void SetInfoField(const std::string& field, const std::vector<s32>& values);
  void SetInfoField(const std::string& field, const std::vector<std::string>& values);
  void AddInfoFieldFlag(const std::string& field);
  void RemoveInfoFieldFlag(const std::string& field);

  template <typename T>
  std::vector<T> GetFormatField(const std::string& field) const;
  void SetFormatField(const std::string& field, const std::vector<f32>& values);
  void SetFormatField(const std::string& field, const std::vector<s32>& values);
  void SetFormatField(const std::string& field, const std::vector<const char*>& values);
  void SetFormatField(const std::string& field, const std::vector<std::string>& values);
  std::string GetGTField() const;
  std::string GetGTField(s32 sample_index) const;
  void SetGTField(const std::string& value);
  void SetGTField(const std::string& value, s32 sample_index);

  void SetChromosome(const std::string& chromosome);
  void SetPosition(s32 position);
  void SetId(const std::string& id);
  void SetAlleles(const std::vector<std::string>& alleles);
  std::optional<f32> GetQuality() const;
  f32 GetQuality(f32 default_value) const;
  void SetQuality(const std::optional<f32>& quality);
  void SetFilter(const char* filter);
  void SetFilter(const std::string& filter);
  void SetFilter(std::string_view filter);
  void AddFilter(const char* filter);
  void AddFilter(const std::string& filter);
  void AddFilter(std::string_view filter);
  void ClearFilters();
  std::vector<std::string> GetFilters() const;
  VcfRecordPtr Clone(const VcfHeaderPtr& new_header) const;

  // used to expose _record to VcfWriter
  friend class VcfWriter;

 private:
  // variables
  BcfHeaderPtr _hdr;
  BcfRecordPtr _record;
  static constexpr char kPhasedAlleleSeparator = '|';
  static constexpr char kUnphasedAlleleSeparator = '/';
  static constexpr char kMissingAllele = '.';
  // Matches a VCF GT string: a single allele (`.` or integer), optionally followed by more
  // alleles joined by `/` (unphased) or `|` (phased), but not a mix of both separators.
  static constexpr auto kGenotypeRegex = R"(^(?:\.|\d+)(?:(?:/(?:\.|\d+))*|(?:\|(?:\.|\d+))*)$)";
};
}  // namespace xoos::io
