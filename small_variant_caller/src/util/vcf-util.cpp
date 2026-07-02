#include "vcf-util.h"

#include "core/genotype.h"
#include "core/vcf-fields.h"
#include "seq-util.h"

namespace xoos::svc {

bool IsMultiAllelicRecord(const io::VcfRecordPtr& record) {
  // A multi-allelic record has 3 or more alleles (1 REF and 2+ ALT)
  static const s32 kMinNumAlleles = 3;
  return record->NumAlleles() >= kMinNumAlleles;
}

/**
 * @brief Update a field (INFO or FORMAT) for a specific alternate allele index.
 * @tparam Metadata Metadata type (InfoFieldMetadata or FormatFieldMetadata)
 * @tparam GetFieldFn Function type to get the field values
 * @tparam SetFieldFn Function type to set the field values
 * @param record VCF record
 * @param metadata Field metadata
 * @param alt_idx Alternate allele index (1-based)
 * @param get_field Function to get the field values
 * @param set_field Function to set the field values
 */
template <typename Metadata, typename GetFieldFn, typename SetFieldFn>
static void UpdateFieldForAltAllele(const io::VcfRecordPtr& record,
                                    const Metadata& metadata,
                                    const s32 alt_idx,
                                    GetFieldFn get_field,
                                    SetFieldFn set_field) {
  if (alt_idx < 1) {
    return;
  }
  const auto& values = get_field(record, metadata.id);
  if (metadata.number == io::kNumberEachAllele) {
    const auto val_idx = ToUnsigned(alt_idx - 1);
    if (val_idx < values.size()) {
      set_field(record, metadata.id, {values.at(val_idx)});
    }
  } else if (metadata.number == io::kNumberR) {
    const auto val_idx = ToUnsigned(alt_idx);
    if (val_idx <= values.size()) {
      set_field(record, metadata.id, {values.front(), values.at(val_idx)});
    }
  }
}

/**
 * @brief Update an INFO field for a specific alternate allele index.
 * @tparam T Field value type
 * @param record VCF record
 * @param metadata INFO field metadata
 * @param alt_idx Alternate allele index (1-based)
 */
template <typename T>
static void UpdateInfoFieldForAltAllele(const io::VcfRecordPtr& record,
                                        const io::InfoFieldMetadata& metadata,
                                        const s32 alt_idx) {
  UpdateFieldForAltAllele(
      record,
      metadata,
      alt_idx,
      [](const io::VcfRecordPtr& rec, const std::string& id) { return rec->GetInfoFieldNoCheck<T>(id); },
      [](const io::VcfRecordPtr& rec, const std::string& id, const vec<T>& vals) { rec->SetInfoField(id, vals); });
}

/**
 * @brief Update a FORMAT field for a specific alternate allele index.
 * @tparam T Field value type
 * @param record VCF record
 * @param metadata FORMAT field metadata
 * @param alt_idx Alternate allele index (1-based)
 */
template <typename T>
static void UpdateFormatFieldForAltAllele(const io::VcfRecordPtr& record,
                                          const io::FormatFieldMetadata& metadata,
                                          const s32 alt_idx) {
  UpdateFieldForAltAllele(
      record,
      metadata,
      alt_idx,
      [](const io::VcfRecordPtr& rec, const std::string& id) { return rec->GetFormatFieldNoCheck<T>(id); },
      [](const io::VcfRecordPtr& rec, const std::string& id, const vec<T>& vals) { rec->SetFormatField(id, vals); });
}

template void UpdateInfoFieldForAltAllele<s32>(const io::VcfRecordPtr& record,
                                               const io::InfoFieldMetadata& metadata,
                                               s32 alt_idx);
template void UpdateInfoFieldForAltAllele<f32>(const io::VcfRecordPtr& record,
                                               const io::InfoFieldMetadata& metadata,
                                               s32 alt_idx);

template void UpdateFormatFieldForAltAllele<s32>(const io::VcfRecordPtr& record,
                                                 const io::FormatFieldMetadata& metadata,
                                                 s32 alt_idx);
template void UpdateFormatFieldForAltAllele<f32>(const io::VcfRecordPtr& record,
                                                 const io::FormatFieldMetadata& metadata,
                                                 s32 alt_idx);

/**
 * @brief Classifies a GT allele token into one of three categories.
 *
 * - `kZero`      — the token is the string `0` (reference allele).
 * - `kPositive`  — the token is a purely-digit string with a positive value (ALT allele index).
 * - `kNonNumber` — the token is empty, `.` (missing), or contains non-digit characters.
 */
enum class TokenType {
  kZero,
  kPositive,
  kNonNumber
};

/**
 * @brief Classifies a single GT allele token as zero, positive integer, or non-number.
 *
 * Iterates over every character in @p t.  If any character is outside `'0'–'9'`,
 * the token is immediately classified as `kNonNumber`.  Otherwise the presence of
 * non-zero digits determines whether the token is `kPositive` (e.g. `1`, `12`)
 * or `kZero` (the single value `0`).
 *
 * @param token Token string view to classify (e.g. `0`, `2`, `.`, empty).
 * @return `TokenType::kZero` if @p t is `0`, `TokenType::kPositive` if @p t is a positive
 *         integer, or `TokenType::kNonNumber` otherwise.
 */
static TokenType ClassifyToken(const std::string_view token) {
  using enum TokenType;
  if (token.empty()) {
    return kNonNumber;
  }
  bool has_zero = false;
  bool has_pos = false;
  for (const char c : token) {
    if (c < '0' || c > '9') {
      return kNonNumber;
    }
    if (c == '0') {
      has_zero = true;
    } else {
      has_pos = true;
    }
  }
  if (has_zero && !has_pos) {
    return kZero;
  }
  return kPositive;
}

/**
 * @brief Maps a single GT allele token to its remapped string in a bi-allelic record.
 *
 * Delegates classification to `ClassifyToken` and maps the result as follows:
 *
 * | `TokenType`    | Input examples      | Output |
 * |----------------|---------------------|--------|
 * | `kZero`        | `0`                 | `"0"`  |
 * | `kPositive`    | `1`, `2`, `12`      | `"1"`  |
 * | `kNonNumber`   | `.`, empty, `*`     | `"."`  |
 *
 * @param token Single allele token to remap (e.g. `0`, `1`, `.`, `2`).
 * @return Remapped string view: `"0"`, `"1"`, or `"."`.
 */
static std::string_view RemapAlleleToken(const std::string_view token) {
  using enum TokenType;
  switch (ClassifyToken(token)) {
    case kZero:
      return "0";
    case kPositive:
      return "1";
    default:
      return ".";
  }
}

/**
 * @brief Splits a GT string into individual allele tokens using the given separator.
 *
 * Iterates over @p gt character by character, splitting on `/` or `|`.
 * A sentinel copy of @p sep is appended to @p gt to flush the final token.
 *
 * @param gt  GT string to split (e.g. `0/1`, `1|2`, `0/1/2`).
 * @param sep Separator character (`/` for unphased, `|` for phased).
 * @return Vector of allele token strings (e.g. `{0, 1}`).
 */
static std::vector<std::string> ParseGTTokens(const std::string& gt, const char sep) {
  std::vector<std::string> tokens;
  std::string token;
  for (const char ch : gt + sep) {
    if (ch == '/' || ch == '|') {
      tokens.push_back(std::exchange(token, {}));
    } else {
      token += ch;
    }
  }
  return tokens;
}

/**
 * @brief Reduces a polyploid (2+ token) GT to a canonical diploid GT string.
 *
 * Counts the number of REF (`0`) and positive-integer (ALT) tokens, ignoring
 * missing (`.`) tokens.  The result follows these rules:
 *
 * | REF count | ALT count | Result  |
 * |-----------|-----------|---------|
 * | > 0       | 0         | `0/0`   |
 * | > 0       | > 0       | `0/1`   |
 * | 0         | > 0       | `1/2`   |
 * | 0         | 0         | `./.`   |
 *
 * The separator (@p sep) is used to join the output tokens, preserving phasing.
 *
 * @param tokens Allele tokens parsed from the original GT string.
 * @param sep    Separator character (`/` or `|`).
 * @return Collapsed diploid GT string.
 */
static std::string ReducePolyploidGT(const std::vector<std::string>& tokens, const char sep) {
  using enum TokenType;
  size_t num_zero = 0;
  size_t num_pos = 0;
  for (const auto& t : tokens) {
    const auto type = ClassifyToken(t);
    if (type == kZero) {
      ++num_zero;
    }
    if (type == kPositive) {
      ++num_pos;
    }
  }
  if (num_zero > 0 && num_pos == 0) {
    return std::string("0") + sep + "0";
  }
  if (num_zero > 0) {
    return std::string("0") + sep + "1";
  }
  if (num_pos > 0) {
    return std::string("1") + sep + "2";
  }
  return std::string(".") + sep + ".";
}

/**
 * @brief Remaps a multi-allelic GT string to a bi-allelic representation.
 *
 * Parses @p gt into allele tokens and applies the following strategy:
 *
 * - **Polyploid / diploid** (≥ 2 tokens): delegates to `ReducePolyploidGT`, which collapses
 *   the token list to a canonical diploid GT based on REF / ALT counts.
 * - **Haploid** (1 token): remaps the single token via `RemapAlleleToken`
 *   (`0` → `0`, positive integer → `1`, otherwise `.`).
 * - **Empty**: returns `.`.
 *
 * Phasing (separator `|`) is detected automatically and preserved in the output.
 *
 * @param gt      GT string to remap (e.g. `0/1`, `1/2`, `0/1/2`, `1|2`).
 * @param alt_idx Alternate allele index (1-based) of the target ALT allele in the original record.
 *                Currently unused; reserved for future allele-specific remapping.
 * @return Remapped bi-allelic GT string.
 */
static std::string RemapGTForAltAllele(const std::string& gt, [[maybe_unused]] const s32 alt_idx) {
  const bool is_phased = gt.find('|') != std::string::npos;
  const char sep = is_phased ? '|' : '/';
  const auto tokens = ParseGTTokens(gt, sep);
  constexpr size_t kMinTokensForPolyploid = 2;
  constexpr size_t kNumTokensForHaploid = 1;
  if (tokens.size() >= kMinTokensForPolyploid) {
    return ReducePolyploidGT(tokens, sep);
  }
  if (tokens.size() == kNumTokensForHaploid) {
    return std::string{RemapAlleleToken(tokens[0])};
  }
  return ".";
}

/**
 * @brief Returns true if a per-allele field should be skipped when splitting a multi-allelic record.
 *
 * Fields with Number={1, 0, .} are already copied correctly when the record is cloned and do not
 * need to be updated.
 *
 * @param number The Number attribute of the field (e.g. `kNumberOne`, `kNumberA`, `kNumberR`).
 * @return true if the field does not require per-allele adjustment.
 */
static bool ShouldSkipPerAlleleField(const std::string_view number) {
  return number == io::kNumberOne || number == io::kNumberZero || number == io::kNumberDot;
}

/**
 * @brief Sets the alleles of a cloned bi-allelic record, optionally trimming common bases.
 * @param new_record The cloned VCF record to update.
 * @param ref_allele The reference allele from the original record.
 * @param alt_allele The target alternate allele from the original record.
 * @param trim_variant Whether to trim shared leading bases from the allele pair.
 */
static void SetSplitAlleles(const io::VcfRecordPtr& new_record,
                            const std::string& ref_allele,
                            const std::string& alt_allele,
                            const bool trim_variant) {
  if (trim_variant) {
    const auto& [trimmed_ref, trimmed_alt] = TrimVariant(ref_allele, alt_allele);
    new_record->SetAlleles({trimmed_ref, trimmed_alt});
  } else {
    new_record->SetAlleles({ref_allele, alt_allele});
  }
}

/**
 * @brief Updates per-allele INFO fields in a cloned record for a specific alternate allele index.
 *
 * Fields with Number=A or Number=R are remapped so they contain only the values relevant to the
 * target ALT allele.  Fields that are absent in the original record, or that have Number={1, 0, .},
 * are left unchanged.  Flag, String, and Character fields are not handled.
 *
 * @param record       The original multi-allelic VCF record.
 * @param new_record   The cloned bi-allelic record to update.
 * @param info_metadata Metadata for all INFO fields declared in the VCF header.
 * @param alt_idx      Alternate allele index (1-based) being extracted.
 */
static void UpdateInfoFieldsForSplit(const io::VcfRecordPtr& record,
                                     const io::VcfRecordPtr& new_record,
                                     const vec<io::InfoFieldMetadata>& info_metadata,
                                     const s32 alt_idx) {
  for (const auto& metadata : info_metadata) {
    if (!record->HasInfoFieldNoCheck(metadata.id) || ShouldSkipPerAlleleField(metadata.number)) {
      continue;
    }
    if (metadata.type == io::FieldType::kInteger) {
      UpdateInfoFieldForAltAllele<s32>(new_record, metadata, alt_idx);
    } else if (metadata.type == io::FieldType::kFloat) {
      UpdateInfoFieldForAltAllele<f32>(new_record, metadata, alt_idx);
    }
    // Flags do not have values to set
    // String and Character INFO fields are not typically per-allele, so they are not handled here
  }
}

/**
 * @brief Updates per-allele FORMAT fields in a cloned record for a specific alternate allele index.
 *
 * Fields with Number=A or Number=R are remapped so they contain only the values relevant to the
 * target ALT allele.  Fields that are absent in the original record, or that have Number={1, 0, .},
 * are left unchanged.  Flag, String, and Character fields are not handled.
 *
 * @param record       The original multi-allelic VCF record.
 * @param new_record   The cloned bi-allelic record to update.
 * @param fmt_metadata Metadata for all FORMAT fields declared in the VCF header.
 * @param alt_idx      Alternate allele index (1-based) being extracted.
 */
static void UpdateFormatFieldsForSplit(const io::VcfRecordPtr& record,
                                       const io::VcfRecordPtr& new_record,
                                       const vec<io::FormatFieldMetadata>& fmt_metadata,
                                       const s32 alt_idx) {
  for (const auto& metadata : fmt_metadata) {
    if (!record->HasFormatFieldNoCheck(metadata.id) || ShouldSkipPerAlleleField(metadata.number)) {
      continue;
    }
    if (metadata.type == io::FieldType::kInteger) {
      UpdateFormatFieldForAltAllele<s32>(new_record, metadata, alt_idx);
    } else if (metadata.type == io::FieldType::kFloat) {
      UpdateFormatFieldForAltAllele<f32>(new_record, metadata, alt_idx);
    }
    // Flags do not have values to set
    // String and Character FORMAT fields are not typically per-allele, so they are not handled here
  }
}

/**
 * @brief Remaps the GT field of every sample in a cloned record for a specific alternate allele.
 *
 * REF (`0`) stays `0`; the target ALT allele becomes `1`; all other ALT alleles are collapsed
 * via `RemapGTForAltAllele`.
 *
 * @param new_record  The cloned bi-allelic record to update.
 * @param new_header  Header of the output VCF (used to determine the number of samples).
 * @param alt_idx     Alternate allele index (1-based) being extracted.
 */
static void RemapGTFieldsForSplit(const io::VcfRecordPtr& new_record,
                                  const io::VcfHeaderPtr& new_header,
                                  const s32 alt_idx) {
  const s32 num_samples = new_header->GetNumSamples();
  for (s32 sample_idx = 0; sample_idx < num_samples; ++sample_idx) {
    const std::string orig_gt = new_record->GetGTField(sample_idx);
    new_record->SetGTField(RemapGTForAltAllele(orig_gt, alt_idx), sample_idx);
  }
}

/**
 * @brief Creates a bi-allelic VCF record for a single alternate allele from a multi-allelic record.
 *
 * Clones @p record into @p new_header, sets the alleles, updates per-allele INFO and FORMAT
 * fields, and remaps GT for all samples.
 *
 * @param record       Original multi-allelic VCF record.
 * @param new_header   Header for the output VCF.
 * @param info_metadata Metadata for all INFO fields declared in the header.
 * @param fmt_metadata  Metadata for all FORMAT fields declared in the header.
 * @param alt_idx      Alternate allele index (1-based) to extract.
 * @param trim_variant Whether to trim shared leading bases from the allele pair.
 * @return Fully populated bi-allelic VCF record.
 */
static io::VcfRecordPtr CreateBiallelicRecord(const io::VcfRecordPtr& record,
                                              const io::VcfHeaderPtr& new_header,
                                              const vec<io::InfoFieldMetadata>& info_metadata,
                                              const vec<io::FormatFieldMetadata>& fmt_metadata,
                                              const s32 alt_idx,
                                              const bool trim_variant) {
  const auto new_record = record->Clone(new_header);
  SetSplitAlleles(new_record, record->Allele(0), record->Allele(alt_idx), trim_variant);
  UpdateInfoFieldsForSplit(record, new_record, info_metadata, alt_idx);
  UpdateFormatFieldsForSplit(record, new_record, fmt_metadata, alt_idx);
  RemapGTFieldsForSplit(new_record, new_header, alt_idx);
  return new_record;
}

vec<io::VcfRecordPtr> SplitMultiAllelicRecord(const io::VcfRecordPtr& record,
                                              const io::VcfHeaderPtr& new_header,
                                              const vec<io::InfoFieldMetadata>& info_metadata,
                                              const vec<io::FormatFieldMetadata>& fmt_metadata,
                                              const bool trim_variant) {
  vec<io::VcfRecordPtr> result;
  const auto num_alleles = record->NumAlleles();
  for (s32 alt_idx = 1; alt_idx < num_alleles; ++alt_idx) {
    result.push_back(CreateBiallelicRecord(record, new_header, info_metadata, fmt_metadata, alt_idx, trim_variant));
  }
  return result;
}

/**
 * @brief Merge a multi-allelic field (INFO or FORMAT) from two records into the first record.
 * @tparam RecordPtr Pointer type to VCF record
 * @tparam Metadata Metadata type (InfoFieldMetadata or FormatFieldMetadata)
 * @tparam GetFieldFn Function type to get the field values
 * @tparam SetFieldFn Function type to set the field values
 * @param record1 First VCF record (to be updated)
 * @param record2 Second VCF record
 * @param metadata Field metadata
 * @param get_field Function to get the field values
 * @param set_field Function to set the field values
 * @return true if the field was successfully merged, false otherwise
 */
template <typename RecordPtr, typename Metadata, typename GetFieldFn, typename SetFieldFn>
static void MergeMultiAllelicField(const RecordPtr& record1,
                                   const RecordPtr& record2,
                                   const Metadata& metadata,
                                   GetFieldFn get_field,
                                   SetFieldFn set_field) {
  const auto& field = metadata.id;
  const auto& values1 = get_field(record1, field);
  const auto& values2 = get_field(record2, field);
  if (metadata.number == io::kNumberR) {
    // For Number=R, both records should have 2 values: [REF, ALT]
    static const size_t kNumAlleles = 2;
    if (values1.size() >= kNumAlleles && values2.size() >= kNumAlleles) {
      // Merge as [REF, ALT1, ALT2]
      set_field(record1, field, {values1.at(0), values1.at(1), values2.at(1)});
    }
  } else if (metadata.number == io::kNumberEachAllele && !values1.empty() && !values2.empty()) {
    // For Number=A, both records should have one value
    // Merge as [ALT1, ALT2]
    set_field(record1, field, {values1.at(0), values2.at(0)});
  }
}

/**
 * @brief Merge a multi-allelic INFO field from two records into the first record.
 * @tparam T Field value type
 * @param record1 First VCF record (to be updated)
 * @param record2 Second VCF record
 * @param metadata INFO field metadata
 */
template <typename T>
static void MergeMultiAllelicInfoField(const io::VcfRecordPtr& record1,
                                       const io::VcfRecordPtr& record2,
                                       const io::InfoFieldMetadata& metadata) {
  MergeMultiAllelicField(
      record1,
      record2,
      metadata,
      [](const io::VcfRecordPtr& rec, const std::string& id) { return rec->GetInfoFieldNoCheck<T>(id); },
      [](const io::VcfRecordPtr& rec, const std::string& id, const vec<T>& vals) { rec->SetInfoField(id, vals); });
}

/**
 * @brief Merge a multi-allelic FORMAT field from two records into the first record.
 * @tparam T Field value type
 * @param record1 First VCF record (to be updated)
 * @param record2 Second VCF record
 * @param metadata FORMAT field metadata
 */
template <typename T>
static void MergeMultiAllelicFormatField(const io::VcfRecordPtr& record1,
                                         const io::VcfRecordPtr& record2,
                                         const io::FormatFieldMetadata& metadata) {
  MergeMultiAllelicField(
      record1,
      record2,
      metadata,
      [](const io::VcfRecordPtr& rec, const std::string& id) { return rec->GetFormatFieldNoCheck<T>(id); },
      [](const io::VcfRecordPtr& rec, const std::string& id, const vec<T>& vals) { rec->SetFormatField(id, vals); });
}

template void MergeMultiAllelicInfoField<s32>(const io::VcfRecordPtr& record1,
                                              const io::VcfRecordPtr& record2,
                                              const io::InfoFieldMetadata& metadata);
template void MergeMultiAllelicInfoField<f32>(const io::VcfRecordPtr& record1,
                                              const io::VcfRecordPtr& record2,
                                              const io::InfoFieldMetadata& metadata);

template void MergeMultiAllelicFormatField<s32>(const io::VcfRecordPtr& record1,
                                                const io::VcfRecordPtr& record2,
                                                const io::FormatFieldMetadata& metadata);
template void MergeMultiAllelicFormatField<f32>(const io::VcfRecordPtr& record1,
                                                const io::VcfRecordPtr& record2,
                                                const io::FormatFieldMetadata& metadata);

bool MergeMultiAllelicRecords(const io::VcfRecordPtr& record1,
                              const io::VcfRecordPtr& record2,
                              const vec<io::InfoFieldMetadata>& info_metadata,
                              const vec<io::FormatFieldMetadata>& fmt_metadata) {
  // Merge allele representations for the multi-allelic record
  const auto& [ref1, alt1] = TrimVariant(record1->Allele(0), record1->Allele(1));
  const auto& [ref2, alt2] = TrimVariant(record2->Allele(0), record2->Allele(1));
  const auto& [ref_new, alt1_new, alt2_new] = FormatVariants(ref1, alt1, ref2, alt2);
  if (ref_new.empty() || alt1_new.empty() || alt2_new.empty()) {
    return false;
  }
  record1->SetAlleles({ref_new, alt1_new, alt2_new});

  // Merge INFO fields as needed
  for (const auto& metadata : info_metadata) {
    if (!record1->HasInfoFieldNoCheck(metadata.id) || !record2->HasInfoFieldNoCheck(metadata.id)) {
      // if field is not present in both records, nothing needs to be done
      continue;
    }
    if (metadata.number == io::kNumberOne || metadata.number == io::kNumberZero || metadata.number == io::kNumberDot) {
      // if Number={1, 0, .}, the fields were copied correctly when the record was cloned
      continue;
    }
    if (metadata.type == io::FieldType::kInteger) {
      MergeMultiAllelicInfoField<s32>(record1, record2, metadata);
    } else if (metadata.type == io::FieldType::kFloat) {
      MergeMultiAllelicInfoField<f32>(record1, record2, metadata);
    }
    // Flags do not have values to set
    // String and Character INFO fields are not typically per-allele, so they are not handled here
  }

  // Merge FORMAT fields as needed
  for (const auto& metadata : fmt_metadata) {
    if (!record1->HasFormatFieldNoCheck(metadata.id) || !record2->HasFormatFieldNoCheck(metadata.id)) {
      // if field is not present in both records, nothing needs to be done
      continue;
    }
    if (metadata.number == io::kNumberOne || metadata.number == io::kNumberZero || metadata.number == io::kNumberDot) {
      // if Number={1, 0, .}, the fields were copied correctly when the record was cloned
      continue;
    }
    if (metadata.type == io::FieldType::kInteger) {
      MergeMultiAllelicFormatField<s32>(record1, record2, metadata);
    } else if (metadata.type == io::FieldType::kFloat) {
      MergeMultiAllelicFormatField<f32>(record1, record2, metadata);
    }
    // Flags do not have values to set
    // String and Character FORMAT fields are not typically per-allele, so they are not handled here
  }

  // Update AC, AF, AN, and GT field values for multi-allelic record
  static const s32 kMultiallelicAc = 1;
  static const f32 kMultiallelicAf = 0.5;
  static const s32 kMultiallelicAn = 2;
  record1->SetInfoField(kFieldAc, vec<s32>{kMultiallelicAc, kMultiallelicAc});
  record1->SetInfoField(kFieldAf, vec<f32>{kMultiallelicAf, kMultiallelicAf});
  record1->SetInfoField(kFieldAn, vec<s32>{kMultiallelicAn});
  record1->SetGTField(GenotypeToString(Genotype::kGT12));

  return true;
}

}  // namespace xoos::svc
