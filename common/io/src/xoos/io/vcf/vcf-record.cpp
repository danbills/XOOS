#include "xoos/io/vcf/vcf-record.h"

#include <numeric>
#include <regex>

#include <htslib/vcf.h>

#include <xoos/error/error.h>
#include <xoos/io/malloc-ptr.h>

#include "xoos/io/vcf/constants.h"

namespace xoos::io {

VcfRecord::VcfRecord(const VcfHeaderPtr& hdr) : _hdr(hdr->_hdr) {
}

/**
 * Create a new record from the header at the specified unpack level.
 * @param hdr Header
 * @param unpack Unpack level
 * @return New record
 */
VcfRecordPtr VcfRecord::CreateFromHeader(const VcfHeaderPtr& hdr, const s32 unpack) {
  auto vcf_record = std::make_shared<VcfRecord>(hdr);
  const auto bcf_record = BcfRecordPtr(bcf_init(), bcf_destroy);
  if (bcf_record == nullptr) {
    throw error::Error("Failed to initialize record");
  }
  vcf_record->_record = bcf_record;
  bcf_unpack(vcf_record->_record.get(), unpack);
  return vcf_record;
}

/**
 * Create a new record from the header at the maximum unpack level.
 * @param hdr Header
 * @return New record
 */
VcfRecordPtr VcfRecord::CreateFromHeader(const VcfHeaderPtr& hdr) {
  return CreateFromHeader(hdr, BCF_UN_ALL);
}

/**
 * Read record from file at the specified unpack level.
 * @param hdr Header
 * @param input_vcf_name HTS file handle
 * @param unpack Unpack level
 * @return Next record
 */
VcfRecordPtr VcfRecord::ReadFromFile(const VcfHeaderPtr& hdr,
                                     const HtsFileSharedPtr& input_vcf_name,
                                     const s32 unpack) {
  const auto record = BcfRecordPtr{bcf_init(), bcf_destroy};
  if (record == nullptr) {
    throw error::Error("Failed to initialize record");
  }
  const auto result = bcf_read(input_vcf_name.get(), hdr->_hdr.get(), record.get());
  if (result == -1) {
    // EOF (normal end of file)
    return nullptr;
  } else if (result < -1) {
    // Report parsing error
    throw error::Error("VCF parsing failed");
  }
  auto shared_record = std::make_shared<VcfRecord>(hdr);
  shared_record->_record = record;
  bcf_unpack(record.get(), unpack);
  return shared_record;
}

/**
 * Read entire record from file.
 * @param hdr Header
 * @param input_vcf_name HTS file handle
 * @return Next record
 */
VcfRecordPtr VcfRecord::ReadFromFile(const VcfHeaderPtr& hdr, const HtsFileSharedPtr& input_vcf_name) {
  return ReadFromFile(hdr, input_vcf_name, BCF_UN_ALL);
}

/**
 * @brief Read the next record in the target region.
 * @param hdr Header
 * @param input_vcf_fp HTS file handle
 * @param tbx_idx Tabix index
 * @param hts_itr HTS iterator
 * @param kstr kstring struct
 * @param unpack Unpack level
 * @return Next record
 */
VcfRecordPtr VcfRecord::ReadFromRegion(const VcfHeaderPtr& hdr,
                                       const HtsFileSharedPtr& input_vcf_fp,
                                       tbx_t* tbx_idx,
                                       hts_itr_t* hts_itr,
                                       Kstring& kstr,
                                       const s32 unpack) {
  const auto record = BcfRecordPtr{bcf_init(), bcf_destroy};
  if (record == nullptr) {
    throw error::Error("Failed to initialize record");
  }
  if (input_vcf_fp->format.format == bcf) {
    const auto bcf_itr_result = bcf_itr_next(input_vcf_fp.get(), hts_itr, record.get());
    if (bcf_itr_result < 0) {
      return nullptr;
    }
  } else {
    const auto tbx_itr_result = tbx_itr_next(input_vcf_fp.get(), tbx_idx, hts_itr, kstr.Get());
    if (tbx_itr_result < 0) {
      return nullptr;
    }
    if (const auto result = vcf_parse(kstr.Get(), hdr->_hdr.get(), record.get()); result < 0) {
      return nullptr;
    }
    kstr.Clear();
  }
  auto shared_record = std::make_shared<VcfRecord>(hdr);
  shared_record->_record = record;
  bcf_unpack(record.get(), unpack);
  return shared_record;
}

std::string VcfRecord::Chromosome() const {
  return {bcf_hdr_id2name(_hdr.get(), _record->rid)};
}

hts_pos_t VcfRecord::Position() const {
  return _record->pos;
}

bool VcfRecord::IsPass() const {
  return bcf_has_filter(_hdr.get(), _record.get(), const_cast<char*>(kFilterPASS)) == 1 ||
         bcf_has_filter(_hdr.get(), _record.get(), const_cast<char*>(kFilterDOT)) == 1;
}

bool VcfRecord::IsSnp() const {
  return bcf_is_snp(_record.get()) != 0;
}

std::string VcfRecord::Id() const {
  return {_record->d.id};
}

std::string VcfRecord::Allele(const s32 which) const {
  return {_record->d.allele[which]};
}

s32 VcfRecord::NumAlleles() const {
  // TODO: this should be u32, check whether using u32 would have any issues
  return static_cast<s32>(_record->n_allele);
}

void VcfRecord::SetChromosome(const std::string& chromosome) {
  _record->rid = bcf_hdr_name2id(_hdr.get(), chromosome.c_str());
}

void VcfRecord::SetPosition(const s32 position) {
  _record->pos = position;
}

void VcfRecord::SetId(const std::string& id) {
  if (bcf_update_id(_hdr.get(), _record.get(), id.c_str()) < 0) {
    throw error::Error("Failed to update id: {}", id);
  }
}

/**
 * Set the alleles for the record.
 * @param alleles Alleles to set, including the reference allele
 */
void VcfRecord::SetAlleles(const std::vector<std::string>& alleles) {
  std::vector<const char*> tmp;
  tmp.reserve(alleles.size());
  for (const auto& value : alleles) {
    tmp.emplace_back(value.c_str());
  }
  if (bcf_update_alleles(_hdr.get(), _record.get(), &tmp[0], static_cast<s32>(tmp.size())) < 0) {
    throw error::Error("Failed to update alleles");
  }
}

/**
 * @brief Set the Quality score for the record.
 * @param quality Quality score, if std::nullopt, sets the quality to missing.
 */
void VcfRecord::SetQuality(const std::optional<f32>& quality) {
  if (quality) {
    _record->qual = *quality;
  } else {
    bcf_float_set_missing(_record->qual);
  }
}

std::optional<f32> VcfRecord::GetQuality() const {
  return (bcf_float_is_missing(_record->qual) != 0) ? std::nullopt : std::make_optional(_record->qual);
}

f32 VcfRecord::GetQuality(f32 default_value) const {
  return GetQuality().value_or(default_value);
}

/**
 * @brief Check whether an INFO field is in the record, without checking header.
 * @param field Field name
 * @return Whether field is found
 */
bool VcfRecord::HasInfoFieldNoCheck(const std::string& field) const {
  return bcf_get_info_id(_record.get(), bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str())) != nullptr;
}

/**
 * @brief Extract the float values of an INFO field, without checking header.
 * @param field Field name
 * @return Vector of float values
 */
template <>
std::vector<f32> VcfRecord::GetInfoFieldNoCheck(const std::string& field) const {
  // look at record without checking header
  s32 num_entries = 0;
  MallocPtr<f32> buffer = nullptr;
  if (bcf_get_info_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_REAL) < 0) {
    return {};  // field not found
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Extract the int values of an INFO field, without checking header.
 * @param field Field name
 * @return Vector of int values
 */
template <>
std::vector<s32> VcfRecord::GetInfoFieldNoCheck(const std::string& field) const {
  // look at record without checking header
  s32 num_entries = 0;
  MallocPtr<s32> buffer = nullptr;
  if (bcf_get_info_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_INT) < 0) {
    return {};  // field not found
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Extract the value string of an INFO field, without checking header.
 * @param field Field name
 * @return value string
 */
std::string VcfRecord::GetInfoFieldStringNoCheck(const std::string& field) const {
  // look at record without checking header
  s32 num_entries = 0;
  MallocPtr<char> buffer = nullptr;
  if (bcf_get_info_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_STR) < 0) {
    return {};  // field not found
  }
  std::string value(buffer.get(), num_entries);

  // when we use bcf_read to set the bcf1_t record, it seems to add a trailing '\0' to format fields. However, when we
  // explicitly set the format fields, it does not add the trailing '\0'. So this may be a hacky way to check if there
  // is a '\0' at the end of the string and remove it if there is.
  if (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
  return value;
}

/**
 * @brief Check whether a FORMAT field is in the record, without checking header.
 * @param field Field name
 * @return Vector of int values
 */
bool VcfRecord::HasFormatFieldNoCheck(const std::string& field) const {
  return bcf_get_fmt_id(_record.get(), bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str())) != nullptr;
}

/**
 * @brief Extract the float values of a FORMAT field, without checking header.
 * @param field Field name
 * @return Vector of float values
 */
template <>
std::vector<f32> VcfRecord::GetFormatFieldNoCheck(const std::string& field) const {
  // look at record without checking header
  s32 num_entries = 0;
  MallocPtr<f32> buffer = nullptr;
  if (bcf_get_format_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_REAL) < 0) {
    return {};
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Extract the int values of a FORMAT field, without checking header.
 * @param field Field name
 * @return Vector of int values
 */
template <>
std::vector<s32> VcfRecord::GetFormatFieldNoCheck(const std::string& field) const {
  // look at record without checking header
  s32 num_entries = 0;
  MallocPtr<s32> buffer = nullptr;
  if (bcf_get_format_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_INT) < 0) {
    return {};
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Extract the value string of a FORMAT field, without checking header.
 * @param field Field name
 * @return value string
 */
std::string VcfRecord::GetFormatFieldStringNoCheck(const std::string& field) const {
  // look at record without checking header
  if (IsGtStandardFormatField(field)) {
    return {GetGTField()};
  }

  s32 num_entries = 0;
  MallocPtr<char> buffer = nullptr;
  if (bcf_get_format_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_STR) < 0) {
    return {};
  }
  std::string value(buffer.get(), num_entries);

  // when we use bcf_read to set the bcf1_t record, it seems to add a trailing '\0' to format fields. However, when we
  // explicitly set the format fields, it does not add the trailing '\0'. So this may be a hacky way to check if there
  // is a '\0' at the end of the string and remove it if there is.
  if (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
  return value;
}

/**
 * @brief Confirm that a field is in the header.
 * @param field Field name
 * @param field_type Field type, e.g. INFO or FORMAT
 */
void VcfRecord::ConfirmFieldInHeader(const std::string& field, const std::string& field_type) const {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw error::Error("{} field {} not found in header", field_type, field);
  }
}

/**
 * @brief Confirm that an INFO field is in the header.
 * @param field Field name
 */
void VcfRecord::ConfirmInfoFieldInHeader(const std::string& field) const {
  ConfirmFieldInHeader(field, "INFO");
}

/**
 * @brief Confirm that a FORMAT field is in the header.
 * @param field Field name
 */
void VcfRecord::ConfirmFormatFieldInHeader(const std::string& field) const {
  ConfirmFieldInHeader(field, "FORMAT");
}

/**
 * @brief Extract the float values of an INFO field.
 * @param field Field name
 */
template <>
std::vector<f32> VcfRecord::GetInfoField(const std::string& field) const {
  ConfirmInfoFieldInHeader(field);
  s32 num_entries = 0;
  MallocPtr<f32> buffer = nullptr;
  if (bcf_get_info_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_REAL) < 0) {
    throw error::Error("Failed to get INFO field {}", field);
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Extract the int values of an INFO field.
 * @param field Field name
 */
template <>
std::vector<s32> VcfRecord::GetInfoField(const std::string& field) const {
  ConfirmInfoFieldInHeader(field);
  s32 num_entries = 0;
  MallocPtr<s32> buffer = nullptr;
  if (bcf_get_info_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_INT) < 0) {
    throw error::Error("Failed to get INFO field {}", field);
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Extract the string values of an INFO field.
 * @param field Field name
 */
template <>
std::vector<std::string> VcfRecord::GetInfoField(const std::string& field) const {
  ConfirmInfoFieldInHeader(field);
  s32 num_entries = 0;
  MallocPtr<char> buffer = nullptr;
  if (bcf_get_info_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_STR) < 0) {
    throw error::Error("Failed to get INFO field {}", field);
  }
  std::string value(buffer.get(), num_entries);

  // when we use bcf_read to set the bcf1_t record, it seems to add a trailing '\0' to format fields. However, when we
  // explicitly set the format fields, it does not add the trailing '\0'. So this may be a hacky way to check if there
  // is a '\0' at the end of the string and remove it if there is.
  if (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
  return {value};
}

/**
 * @brief Extract the bool value of an INFO field.
 * For compatibility reasons, this will return a vector containing exactly one bool
 * @param field Field name
 */
template <>
std::vector<bool> VcfRecord::GetInfoField(const std::string& field) const {
  ConfirmInfoFieldInHeader(field);
  s32 num_entries = 0;
  // this will be unused for retrieving bool
  MallocPtr<char> buffer = nullptr;
  // from htslib documentation:
  // Returns negative value on error or the number of values (including missing
  // values) put in *dst on success.
  // ...
  // bcf_get_info_flag() does not store anything in *dst but returns 1 if the
  // flag is set or 0 if not.
  return {(bcf_get_info_values(
               _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_FLAG) >
           0)};
}

/**
 * @brief Set the float values of an INFO field.
 * @param field Field name
 * @param values Vector of float values
 */
void VcfRecord::SetInfoField(const std::string& field, const std::vector<f32>& values) {
  ConfirmInfoFieldInHeader(field);
  if (bcf_update_info_float(_hdr.get(), _record.get(), field.c_str(), &values[0], values.size()) < 0) {
    throw error::Error("Failed to update INFO field {}", field);
  }
}

/**
 * @brief Set the int values of an INFO field.
 * @param field Field name
 * @param values Vector of int values
 */
void VcfRecord::SetInfoField(const std::string& field, const std::vector<s32>& values) {
  ConfirmInfoFieldInHeader(field);
  if (bcf_update_info_int32(_hdr.get(), _record.get(), field.c_str(), &values[0], values.size()) < 0) {
    throw error::Error("Failed to update INFO field {}", field);
  }
}

/**
 * @brief Set the string values of an INFO field.
 * @param field Field name
 * @param values Vector of string values
 */
void VcfRecord::SetInfoField(const std::string& field, const std::vector<std::string>& values) {
  ConfirmInfoFieldInHeader(field);
  const auto value = std::accumulate(values.begin(), values.end(), std::string{});
  if (bcf_update_info_string(_hdr.get(), _record.get(), field.c_str(), value.c_str()) < 0) {
    throw error::Error("Failed to update INFO field {}", field);
  }
}

/**
 * @brief Add a flag INFO field to the record.
 * @param field Field name
 */
void VcfRecord::AddInfoFieldFlag(const std::string& field) {
  ConfirmInfoFieldInHeader(field);
  // set to 1 to add this flag
  if (bcf_update_info_flag(_hdr.get(), _record.get(), field.c_str(), nullptr, 1) < 0) {
    throw error::Error("Failed to add INFO field {}", field);
  }
}

/**
 * @brief Remove a flag INFO field from the record.
 * @param field Field name
 */
void VcfRecord::RemoveInfoFieldFlag(const std::string& field) {
  ConfirmInfoFieldInHeader(field);
  // set to 0 to remove this flag
  if (bcf_update_info_flag(_hdr.get(), _record.get(), field.c_str(), nullptr, 0) < 0) {
    throw error::Error("Failed to remove INFO field {}", field);
  }
}

/**
 * @brief Get the float value of an FORMAT field.
 * @param field Field name
 * @return Vector of float values
 */
template <>
std::vector<f32> VcfRecord::GetFormatField(const std::string& field) const {
  ConfirmFormatFieldInHeader(field);
  s32 num_entries = 0;
  MallocPtr<f32> buffer = nullptr;
  if (bcf_get_format_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_REAL) < 0) {
    throw error::Error("Failed to get FORMAT field {}", field);
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Get the int value of an FORMAT field.
 * @param field Field name
 * @return Vector of int values
 */
template <>
std::vector<s32> VcfRecord::GetFormatField(const std::string& field) const {
  ConfirmFormatFieldInHeader(field);
  s32 num_entries = 0;
  MallocPtr<s32> buffer = nullptr;
  if (bcf_get_format_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_INT) < 0) {
    throw error::Error("Failed to get FORMAT field {}", field);
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Get the string value of an FORMAT field.
 * @param field Field name
 * @return Vector of string values
 */
template <>
std::vector<std::string> VcfRecord::GetFormatField(const std::string& field) const {
  ConfirmFormatFieldInHeader(field);

  if (IsGtStandardFormatField(field)) {
    return {GetGTField()};
  }

  s32 num_entries = 0;
  MallocPtr<char> buffer = nullptr;
  if (bcf_get_format_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_STR) < 0) {
    throw error::Error("Failed to get FORMAT field {}", field);
  }
  std::string value(buffer.get(), num_entries);

  // when we use bcf_read to set the bcf1_t record, it seems to add a trailing '\0' to format fields. However, when we
  // explicitly set the format fields, it does not add the trailing '\0'. So this may be a hacky way to check if there
  // is a '\0' at the end of the string and remove it if there is.
  if (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
  return {value};
}

/**
 * @brief Set the float values of a FORMAT field.
 * @param field Field name
 * @param values Vector of float values
 */
void VcfRecord::SetFormatField(const std::string& field, const std::vector<f32>& values) {
  ConfirmFormatFieldInHeader(field);
  if (bcf_update_format_float(_hdr.get(), _record.get(), field.c_str(), &values[0], values.size()) < 0) {
    throw error::Error("Failed to update FORMAT field {}", field);
  }
}

/**
 * @brief Set the int values of a FORMAT field.
 * @param field Field name
 * @param values Vector of int values
 */
void VcfRecord::SetFormatField(const std::string& field, const std::vector<s32>& values) {
  ConfirmFormatFieldInHeader(field);
  if (bcf_update_format_int32(_hdr.get(), _record.get(), field.c_str(), &values[0], values.size()) < 0) {
    throw error::Error("Failed to update FORMAT field {}", field);
  }
}

/**
 * @brief Set the string values of a FORMAT field.
 * @param field Field name
 * @param values Vector of C-string values
 * @note If the field is GT, it will be set using SetGTField instead.
 */
void VcfRecord::SetFormatField(const std::string& field, const std::vector<const char*>& values) {
  ConfirmFormatFieldInHeader(field);

  if (IsGtStandardFormatField(field)) {
    SetGTField(values[0]);
    return;
  }

  std::vector tmp(values);
  if (bcf_update_format_string(_hdr.get(), _record.get(), field.c_str(), &tmp[0], static_cast<s32>(tmp.size())) < 0) {
    throw error::Error("Failed to update FORMAT field {}", field);
  }
}

/**
 * @brief Set the string values of a FORMAT field.
 * @param field Field name
 * @param values Vector of string values
 * @note If the field is GT, it will be set using SetGTField instead.
 */
void VcfRecord::SetFormatField(const std::string& field, const std::vector<std::string>& values) {
  ConfirmFormatFieldInHeader(field);

  if (IsGtStandardFormatField(field)) {
    SetGTField(values[0]);
    return;
  }

  std::vector<const char*> tmp;
  tmp.reserve(values.size());
  for (const auto& value : values) {
    tmp.emplace_back(value.c_str());
  }
  if (bcf_update_format_string(_hdr.get(), _record.get(), field.c_str(), &tmp[0], static_cast<s32>(tmp.size())) < 0) {
    throw error::Error("Failed to update FORMAT field {}", field);
  }
}

/**
 * @brief Set the genotype field (GT) in the record for a specific sample.
 * @param value Genotype string, e.g. "0/1" or "1|0"
 * @param sample_index Sample index (0-based)
 * @throws std::runtime_error if the value is invalid or sample_index is out of range
 */
void VcfRecord::SetGTField(const std::string& value, const s32 sample_index) {
  static const std::regex kPattern(kGenotypeRegex);
  if (!std::regex_match(value, kPattern)) {
    throw error::Error("Invalid genotype field value: {}", value);
  }

  // Detect separator: phased (|) or unphased (/)
  const bool is_phased = value.find(kPhasedAlleleSeparator) != std::string::npos;
  const char sep = is_phased ? kPhasedAlleleSeparator : kUnphasedAlleleSeparator;

  // Tokenise on the separator and convert each token to an htslib GT integer
  auto parse_allele = [&is_phased](const std::string& allele_str) {
    if (allele_str == ".") {
      return bcf_gt_missing;
    }
    const s32 allele = std::stoi(allele_str);
    return is_phased ? bcf_gt_phased(allele) : bcf_gt_unphased(allele);
  };

  // Walk the GT string, splitting on the detected separator to build the allele integer vector.
  std::vector<s32> alleles;
  std::string token;
  for (const char ch : value) {
    if (ch == sep) {
      alleles.push_back(parse_allele(token));
      token.clear();
    } else {
      token += ch;
    }
  }
  alleles.push_back(parse_allele(token));
  const auto new_ploidy = alleles.size();

  const s32 num_samples = bcf_hdr_nsamples(_hdr.get());
  if (sample_index < 0 || sample_index >= num_samples) {
    throw error::Error("Sample index out of range: {}", sample_index);
  }
  const auto num_samples_u = static_cast<size_t>(num_samples);
  const auto sample_idx_u = static_cast<size_t>(sample_index);

  // Read the existing GT buffer (needed to preserve other samples' genotypes).
  MallocPtr<s32> existing_gts = nullptr;
  s32 existing_num_gts = 0;
  const bool has_existing =
      bcf_get_genotypes(_hdr.get(), _record.get(), &existing_gts, &existing_num_gts) > 0 && existing_num_gts > 0;

  // Use the larger of the existing ploidy and the new ploidy as the buffer ploidy. This avoids
  // silently truncating other samples' genotype data when the new GT has fewer alleles.
  const size_t existing_ploidy = has_existing ? static_cast<size_t>(existing_num_gts) / num_samples_u : new_ploidy;
  const size_t buf_ploidy = std::max(existing_ploidy, new_ploidy);
  const size_t total_gts = num_samples_u * buf_ploidy;

  // Initialize all positions to vector_end so padding slots beyond each sample's true ploidy
  // are handled correctly by htslib (they serialize as absent rather than as missing alleles).
  auto genotypes = MallocPtr<s32>(static_cast<s32*>(malloc(total_gts * sizeof(s32))));
  std::fill_n(genotypes.get(), total_gts, bcf_int32_vector_end);

  // Copy existing genotypes for all samples except the target, preserving their full allele vectors.
  if (has_existing) {
    for (size_t s = 0; s < num_samples_u; ++s) {
      if (s == sample_idx_u) {
        continue;
      }
      for (size_t a = 0; a < existing_ploidy; ++a) {
        genotypes.get()[s * buf_ploidy + a] = existing_gts.get()[s * existing_ploidy + a];
      }
      // Positions [existing_ploidy, buf_ploidy) remain bcf_int32_vector_end.
    }
  }

  // Write the target sample's new alleles.
  for (size_t i = 0; i < new_ploidy; ++i) {
    genotypes.get()[sample_idx_u * buf_ploidy + i] = alleles[i];
  }
  // Place a vector_end sentinel after the last allele when the new ploidy is shorter than the
  // buffer ploidy, so htslib serializes this sample as e.g. "0/1" instead of "0/1/.".
  if (new_ploidy < buf_ploidy) {
    genotypes.get()[sample_idx_u * buf_ploidy + new_ploidy] = bcf_int32_vector_end;
  }

  if (bcf_update_genotypes(_hdr.get(), _record.get(), genotypes.get(), static_cast<s32>(total_gts)) < 0) {
    throw error::Error("Failed to update FORMAT field GT");
  }
}

/**
 * @brief Set the genotype field (GT) in the record.
 * @param value Genotype string, e.g. "0/1" or "1|0"
 * @throws std::runtime_error if the value is invalid
 */
void VcfRecord::SetGTField(const std::string& value) {
  SetGTField(value, 0);
}

/**
 * @brief Get the genotype field (GT) from the record.
 * @return Genotype string, e.g. "0/1" or "1|0"
 * @throws std::runtime_error if the GT field cannot be retrieved
 */
std::string VcfRecord::GetGTField() const {
  return GetGTField(0);
}

/**
 * @brief Get the genotype field (GT) for a specific sample.
 * @param sample_index Sample index (0-based)
 * @return Genotype string, e.g. "0/1" or "1|0"
 * @throws std::runtime_error if the GT field cannot be retrieved or sample_index is out of range
 */
std::string VcfRecord::GetGTField(const s32 sample_index) const {
  const s32 num_samples = bcf_hdr_nsamples(_hdr.get());
  if (sample_index < 0 || sample_index >= num_samples) {
    throw error::Error("Sample index out of range: {}", sample_index);
  }

  s32 num_entries = 0;
  MallocPtr<s32> buffer = nullptr;
  if (bcf_get_genotypes(_hdr.get(), _record.get(), &buffer, &num_entries) < 0) {
    throw error::Error("Failed to get FORMAT field GT");
  }

  const auto ploidy = static_cast<size_t>(num_entries / num_samples);
  const auto offset = static_cast<size_t>(sample_index) * ploidy;

  // Copy the raw GT integers into a vector to enable bounds-checked access via .at().
  const auto gts = std::vector<s32>(buffer.get(), buffer.get() + static_cast<size_t>(num_entries));

  std::string value;
  for (size_t i = 0; i < ploidy; ++i) {
    if (gts.at(offset + i) == bcf_int32_vector_end) {
      break;
    }
    if (i > 0) {
      value += bcf_gt_is_phased(gts.at(offset + i)) ? kPhasedAlleleSeparator : kUnphasedAlleleSeparator;
    }
    if (bcf_gt_is_missing(gts.at(offset + i))) {
      value += kMissingAllele;
    } else {
      value += std::to_string(bcf_gt_allele(gts.at(offset + i)));
    }
  }
  return value;
}

/**
 * @brief Get the genotype field (GT) for a specific sample.
 * @param filter Filter name
 * @throws std::runtime_error if the filter is not found in the header or if updating fails
 */
void VcfRecord::SetFilter(const char* filter) {
  s32 index = bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, filter);
  if (index < 0) {
    throw error::Error("FILTER {} not found in header", filter);
  }
  if (bcf_update_filter(_hdr.get(), _record.get(), &index, 1) < 0) {
    throw error::Error("Failed to update FILTER {}", filter);
  }
}

/**
 * @brief Set the filter for the record.
 * @param filter Filter name
 * @throws std::runtime_error if the filter is not found in the header or if updating fails
 */
void VcfRecord::SetFilter(const std::string& filter) {
  SetFilter(filter.c_str());
}

/**
 * @brief Set the filter for the record.
 * @param filter Filter name
 * @throws std::runtime_error if the filter is not found in the header or if updating fails
 */
void VcfRecord::SetFilter(const std::string_view filter) {
  const std::string filter_str(filter);
  SetFilter(filter_str.c_str());
}

/**
 * @brief Add a filter to the record.
 * @param filter Filter name
 * @throws std::runtime_error if the filter is not found in the header or if adding fails
 */
void VcfRecord::AddFilter(const char* filter) {
  const s32 index = bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, filter);
  if (index < 0) {
    throw error::Error("FILTER {} not found in header", filter);
  }
  if (bcf_add_filter(_hdr.get(), _record.get(), index) < 0) {
    throw error::Error("Failed to add FILTER {}", filter);
  }
}

/**
 * @brief Add a filter to the record.
 * @param filter Filter name
 * @throws std::runtime_error if the filter is not found in the header or if adding fails
 */
void VcfRecord::AddFilter(const std::string& filter) {
  AddFilter(filter.c_str());
}

/**
 * @brief Add a filter to the record.
 * @param filter Filter name
 * @throws std::runtime_error if the filter is not found in the header or if adding fails
 */
void VcfRecord::AddFilter(const std::string_view filter) {
  const std::string filter_str(filter);
  AddFilter(filter_str.c_str());
}

/**
 * @brief Remove all filters from the record.
 * @throws std::runtime_error if clearing filters fails
 */
void VcfRecord::ClearFilters() {
  if (bcf_update_filter(_hdr.get(), _record.get(), nullptr, 0) < 0) {
    throw error::Error("Failed to clear FILTER");
  }
}

/**
 * @brief Get the filters applied to the record.
 * @return Vector of filter names
 */
std::vector<std::string> VcfRecord::GetFilters() const {
  std::vector<std::string> filters;

  if (_record->d.n_flt > 0) {
    filters.reserve(_record->d.n_flt);
    for (s32 i = 0; i < _record->d.n_flt; ++i) {
      filters.emplace_back(bcf_hdr_int2id(_hdr, BCF_DT_ID, _record->d.flt[i]));
    }
  }
  return filters;
}

/**
 * @brief Clones record with a new header
 * @param new_header
 * @return Cloned record
 */
VcfRecordPtr VcfRecord::Clone(const VcfHeaderPtr& new_header) const {
  const auto record = BcfRecordPtr{bcf_dup(_record.get()), bcf_destroy};
  if (record == nullptr) {
    throw error::Error("Failed to clone record");
  }
  auto shared_record = std::make_shared<VcfRecord>(new_header);
  shared_record->_record = record;
  return shared_record;
}

}  // namespace xoos::io
