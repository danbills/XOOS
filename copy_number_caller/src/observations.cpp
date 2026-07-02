#include "observations.h"

#include <fstream>
#include <set>
#include <sstream>
#include <utility>

#include <csv.hpp>

#include <xoos/error/error.h>
#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "io/column-names.h"
#include "io/write-bigwig.h"
#include "utility/utility-functions.h"

extern "C" {
#include <libbigwig/bigWig.h>
}

namespace xoos::cnc {
/**
 * @brief read a BED file containing observation values (e.g. log ratios, mapping quality).
 * Expected BED format: Contig\tStart\tEnd\tValue (0-based half-open coordinates).
 * `#`-prefixed metadata/comments and a `#`-prefixed column header are supported and skipped.
 * @param in_stream input stream for the BED file
 * @return Observations object
 */
Observations ReadObservations(std::istream& in_stream) {
  Observations ret;
  std::vector<f64> obvs;
  std::vector<arma::uword> starts;
  std::vector<arma::uword> ends;
  std::string line;
  while (std::getline(in_stream, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream line_ss(line);
    std::string contig;
    std::string start_str;
    std::string end_str;
    std::string value_str;
    if (!std::getline(line_ss, contig, '\t') || !std::getline(line_ss, start_str, '\t') ||
        !std::getline(line_ss, end_str, '\t') || !std::getline(line_ss, value_str, '\t')) {
      throw error::Error("BED file must have at least 4 columns (Contig, Start, End, Value): \"{}\"", line);
    }
    arma::uword start, end;
    f64 value;
    try {
      start = std::stoul(start_str);
      end = std::stoul(end_str);
    } catch (const std::exception& e) {
      throw error::Error("BED file parse error: {} (line: \"{}\")", e.what(), line);
    }
    try {
      value = std::stod(value_str);
    } catch (const std::exception&) {
      value = NAN;
    }
    ret.contigs.push_back(contig);
    starts.push_back(start);
    ends.push_back(end);
    obvs.push_back(value);
  }

  ret.obvs = arma::vec(obvs);
  ret.starts = arma::uvec(starts);
  ret.ends = arma::uvec(ends);
  return ret;
}

const std::string kPassFilter = "PASS";

/**
 * @brief check if a VcfRecord has the "PASS" filter
 * @param rec  vcf record
 * @return  bool
 */
static bool VcfRecordIsPass(const io::VcfRecordPtr& rec) {
  std::vector<std::string> filters = rec->GetFilters();
  return std::ranges::find(filters, kPassFilter) != filters.end();
}

/**
 * @brief parse a VCF and return the BAF value for each variant. If the VCF is from SVC, apply additional filters
 * @param vcf_rdr VcfReader object
 * @param keep_all_variants if true, additional filters appropriate for if cell_line mode is enabled
 * @param update_fn function that takes in: VcfRecordPtr, ref_ad, alt_ad, baf
 * @param baf_filter_options BAF filter options
 * @param sample_metadata_options sample metadata options containing the normal sample name used to resolve the
 * correct sample column in potentially multi-sample VCFs
 */
template <typename UpdateFn>
static void AccessBAFValuesFromVCF(io::VcfReader& vcf_rdr,
                                   bool keep_all_variants,
                                   UpdateFn&& update_fn,
                                   const BafFilterOptions& baf_filter_options,
                                   const SampleMetadataOptions& sample_metadata_options) {
  auto hdr = vcf_rdr.GetHeader();
  const auto sample_info = GetSingleSampleInfoFromHeader(hdr, sample_metadata_options.normal_sample_name.value());
  io::VcfRecordPtr rec;
  while ((rec = vcf_rdr.GetNextRecord()) != nullptr) {
    // continue if rec not snv
    if (!rec->IsSnp()) {
      continue;
    }
    // continue if rec has multiple alternative alleles
    if (rec->NumAlleles() > 2) {
      continue;
    }
    // continue if rec does not have PASS filter
    if (!VcfRecordIsPass(rec)) {
      continue;
    }
    // continue if phased
    std::string gt = rec->GetGTField();
    if (gt.find('|') != std::string::npos) {
      continue;
    }
    // TODO: continue if the SNP overlaps a "PASS" deletion
    // get AD field
    auto ad = rec->GetFormatField<s32>("AD");
    if (ad.empty()) {
      Logging::Error("no AD field found!");
      throw std::runtime_error("required VCF field AD missing");
    }

    const auto min_ad_size = 2 * static_cast<size_t>(sample_info.normal_index) + 2;
    if (ad.size() < min_ad_size) {
      Logging::Error("AD field has {} elements but expected at least {} for sample index {}",
                     ad.size(),
                     min_ad_size,
                     sample_info.normal_index);
      throw std::runtime_error("AD field does not contain enough elements for the requested sample");
    }

    // ID column of VCF file may not be set. As such, we will create a mutation ID for each record using the format of
    // CHROM:POS_REF>ALT. This ensures log messages of variants being filtered will be informative.
    auto mutation_id =
        rec->Chromosome() + ":" + std::to_string(rec->Position() + 1) + "_" + rec->Allele(0) + ">" + rec->Allele(1);
    rec->SetId(mutation_id);

    auto ref_ad = ad[2 * sample_info.normal_index];
    auto alt_ad = ad[2 * sample_info.normal_index + 1];
    if (alt_ad + ref_ad == 0) {
      Logging::Debug("rec ID {} has total allele depth of 0. skipping...", rec->Id());
      continue;
    }
    if (alt_ad + ref_ad < baf_filter_options.normal_sample_min_depth) {
      Logging::Debug(
          "rec ID {} has total allele depth < {}. skipping...", rec->Id(), baf_filter_options.normal_sample_min_depth);
      continue;
    }
    f64 baf = static_cast<f64>(alt_ad) / static_cast<f64>(alt_ad + ref_ad);
    if (keep_all_variants ||
        (baf > baf_filter_options.normal_sample_min_baf && baf < baf_filter_options.normal_sample_max_baf)) {
      update_fn(rec, ref_ad, alt_ad, baf);
    } else {
      Logging::Debug("skipping rec {} because BAF is < {} or > {}",
                     rec->Id(),
                     baf_filter_options.normal_sample_min_baf,
                     baf_filter_options.normal_sample_max_baf);
    }
  }
}

/**
 * @brief Parse a VCF and return BAF values for each variant. Assumes the VCF has a tumor and normal column.
 *
 * Supports two VCF types:
 *   - Non-somatic flagged VCF: Unfiltered T-N VCF where the FILTER column is ignored. All biallelic SNPs passing
 * depth/BAF filters are used as germline heterozygous sites. No somatic VAFs are collected.
 *   - Somatic flagged VCF: Unfiltered T-N VCF with flagged somatic variants (FILTER == PASS indicates somatic).
 * Germline het SNPs have a non-PASS filter. The presence of any PASS variant identifies this type. A minimum fraction
 * of heterozygous SNPs is enforced to ensure the VCF contains sufficient germline variants.
 *
 * @param vcf_rdr VcfReader object
 * @param update_fn function that takes in: VcfRecordPtr, ref_ad, alt_ad, baf
 */
template <typename UpdateFn>
static VcfCheckResult AccessBAFValuesFromTumorNormalVCF(io::VcfReader& vcf_rdr,
                                                        UpdateFn&& update_fn,
                                                        const BafFilterOptions& baf_filter_options,
                                                        const SampleMetadataOptions& sample_metadata_options,
                                                        vec<f64>& somatic_vafs) {
  auto sample_info = GetTumorNormalSampleInfoFromHeader(vcf_rdr.GetHeader(),
                                                        sample_metadata_options.normal_sample_name.value(),
                                                        sample_metadata_options.tumor_sample_name.value());
  xoos::io::VcfRecordPtr rec;

  s32 normal_ref_ad;
  s32 normal_alt_ad;
  f64 normal_baf;
  s32 tumor_ref_ad;
  s32 tumor_alt_ad;
  f64 tumor_baf;

  u32 num_somatic_variants = 0;
  u32 num_het_snps = 0;

  while ((rec = vcf_rdr.GetNextRecord()) != nullptr) {
    // continue if rec not snv
    if (!rec->IsSnp()) {
      continue;
    }

    // continue if rec has multiple alternative alleles
    if (rec->NumAlleles() > 2) {
      continue;
    }

    // There are AD values for both REF and ALT
    // REF values are at indexes 0 and 2; ALT values are at indexes 1 and 3
    // throw an error if the rec lacks an FORMAT::AD field
    const auto& ad_vals = rec->GetFormatField<s32>("AD");
    if (ad_vals.empty()) {
      Logging::Error("no AD field found!");
      throw std::runtime_error("required VCF field AD missing");
    }

    const auto min_ad_size = 2 * static_cast<size_t>(std::max(sample_info.normal_index, sample_info.tumor_index)) + 2;
    if (ad_vals.size() < min_ad_size) {
      Logging::Error("AD field has {} elements but expected at least {} for normal index {} and tumor index {}",
                     ad_vals.size(),
                     min_ad_size,
                     sample_info.normal_index,
                     sample_info.tumor_index);
      throw std::runtime_error("AD field does not contain enough elements for the requested samples");
    }

    normal_ref_ad = ad_vals[2 * sample_info.normal_index];      // +0 for REF
    normal_alt_ad = ad_vals[2 * sample_info.normal_index + 1];  // +1 for ALT
    tumor_ref_ad = ad_vals[2 * sample_info.tumor_index];        // +0 for REF
    tumor_alt_ad = ad_vals[2 * sample_info.tumor_index + 1];    // +1 for ALT

    // special handling for somatic flagged variants must happen before filtering
    const auto is_pass = VcfRecordIsPass(rec);
    if (is_pass) {
      // A PASS variant identifies this as a somatic-flagged VCF. Count it regardless of depth
      // so VcfCheckResult::is_sufficient classification is correct.
      ++num_somatic_variants;
      if (tumor_ref_ad + tumor_alt_ad > 0) {
        tumor_baf = static_cast<f64>(tumor_alt_ad) / (tumor_ref_ad + tumor_alt_ad);
        somatic_vafs.push_back(tumor_baf);
      }
      continue;
    }

    // conditions for filtering
    const auto mutation_id =
        rec->Chromosome() + ":" + std::to_string(rec->Position() + 1) + "_" + rec->Allele(0) + ">" + rec->Allele(1);
    const auto normal_total_depth = normal_ref_ad + normal_alt_ad;
    if (normal_total_depth == 0 || normal_total_depth < baf_filter_options.normal_sample_min_depth) {
      Logging::Debug("rec ID {} has total allele depth {} (min {}) in the normal sample {}, skipping",
                     mutation_id,
                     normal_total_depth,
                     baf_filter_options.normal_sample_min_depth,
                     sample_info.normal_sample);
      continue;
    }
    normal_baf = static_cast<f64>(normal_alt_ad) / (normal_ref_ad + normal_alt_ad);
    if (normal_baf < baf_filter_options.normal_sample_min_baf ||
        normal_baf > baf_filter_options.normal_sample_max_baf) {
      Logging::Debug("rec ID {} has BAF < {} or > {} in the normal sample {}",
                     mutation_id,
                     baf_filter_options.normal_sample_min_baf,
                     baf_filter_options.normal_sample_max_baf,
                     sample_info.normal_sample);
      continue;
    }
    if (tumor_ref_ad + tumor_alt_ad == 0) {
      Logging::Debug("rec ID {} has 0 depth in the tumor sample {}", mutation_id, sample_info.tumour_sample);
      continue;
    }
    tumor_baf = static_cast<f64>(tumor_alt_ad) / (tumor_ref_ad + tumor_alt_ad);

    ++num_het_snps;
    update_fn(rec, tumor_ref_ad, tumor_alt_ad, tumor_baf);
  }

  // Build a VcfCheckResult so the caller can decide how to handle data-quality issues.
  const auto num_relevant_variants = num_het_snps + num_somatic_variants;
  const std::optional<f64> het_snp_fraction =
      num_relevant_variants > 0
          ? std::make_optional(static_cast<f64>(num_het_snps) / static_cast<f64>(num_relevant_variants))
          : std::nullopt;
  // The VCF is sufficient when either:
  //   - no somatic variants were detected (not a somatic-flagged VCF), or
  //   - the het-SNP fraction meets the minimum threshold, or
  //   - the user explicitly forced somatic variant parsing.
  const bool is_sufficient = num_somatic_variants == 0 ||
                             (het_snp_fraction.has_value() && *het_snp_fraction >= kDefaultMinHetSnpFraction) ||
                             baf_filter_options.force_enable_somatic_variant_parsing;

  return VcfCheckResult{
      .num_het_snps = num_het_snps,
      .num_somatic_variants = num_somatic_variants,
      .het_snp_fraction = het_snp_fraction,
      .is_sufficient = is_sufficient,
  };
}

// overload for when somatic VAFs are not of interest
template <typename UpdateFn>
static VcfCheckResult AccessBAFValuesFromTumorNormalVCF(io::VcfReader& vcf_rdr,
                                                        UpdateFn&& update_fn,
                                                        const BafFilterOptions& baf_filter_options,
                                                        const SampleMetadataOptions& sample_metadata_options) {
  vec<f64> unused_somatic_vafs;
  return AccessBAFValuesFromTumorNormalVCF(
      vcf_rdr, std::forward<UpdateFn>(update_fn), baf_filter_options, sample_metadata_options, unused_somatic_vafs);
}

/**
 * @brief read a VCF file, get the alt-allele depth  return them in the form of a
 * RefAltObservations object
 * @param fname vcf file name
 * @param keep_all_variants keep homozygous variants
 * @return Observations object
 */
RefAltObservations GetDepthsFromVcf(io::VcfReader& vcf_rdr,
                                    bool keep_all_variants,
                                    bool has_tumor_normal,
                                    const BafFilterOptions& baf_filter_options,
                                    const SampleMetadataOptions& sample_metadata_options) {
  Observations ref_obvs;
  Observations alt_obvs;
  std::vector<arma::uword> starts;
  std::vector<arma::uword> ends;
  std::vector<f64> alt_ads;
  std::vector<f64> ref_ads;
  std::vector<f64> somatic_vafs;
  auto update_fn =
      [&ref_obvs, &alt_obvs, &starts, &ends, &ref_ads, &alt_ads](
          const io::VcfRecordPtr& rec, const s32 ref_ad, const s32 alt_ad, [[maybe_unused]] const f64 baf) {
        ref_obvs.contigs.push_back(rec->Chromosome());
        alt_obvs.contigs.push_back(rec->Chromosome());
        starts.push_back(rec->Position());
        ends.push_back(rec->Position() + 1);
        ref_ads.push_back(static_cast<f64>(ref_ad));
        alt_ads.push_back(static_cast<f64>(alt_ad));
      };
  std::optional<VcfCheckResult> vcf_check;
  if (has_tumor_normal) {
    if (!sample_metadata_options.normal_sample_name.has_value() ||
        !sample_metadata_options.tumor_sample_name.has_value()) {
      throw std::runtime_error("normal_sample_name and tumor_sample_name must be provided if in somatic T/N mode");
    }
    vcf_check = AccessBAFValuesFromTumorNormalVCF(
        vcf_rdr, update_fn, baf_filter_options, sample_metadata_options, somatic_vafs);
  } else {
    if (!sample_metadata_options.normal_sample_name.has_value()) {
      throw std::runtime_error("sample_name must be provided for single-sample VCF parsing");
    }
    AccessBAFValuesFromVCF(vcf_rdr, keep_all_variants, update_fn, baf_filter_options, sample_metadata_options);
  }
  ref_obvs.starts = arma::uvec(starts);
  ref_obvs.ends = arma::uvec(ends);
  ref_obvs.obvs = arma::vec(ref_ads);
  alt_obvs.starts = arma::uvec(starts);
  alt_obvs.ends = arma::uvec(ends);
  alt_obvs.obvs = arma::vec(alt_ads);
  return RefAltObservations{
      ref_obvs, alt_obvs, somatic_vafs.empty() ? std::nullopt : std::make_optional(somatic_vafs), vcf_check};
}

Observations GetDhFromDepths(const Observations& ref_depths, const Observations& alt_depths) {
  if (ref_depths.contigs.size() != alt_depths.contigs.size()) {
    Logging::Error("ref_depths and alt_depths must have the same number of observations");
    throw std::runtime_error("ref_depths and alt_depths must have the same number of observations");
  }
  Observations dh_obvs;
  std::vector<f64> dh_vals;
  std::vector<s32> dps;
  std::vector<arma::uword> keep_starts;
  std::vector<arma::uword> keep_ends;
  for (arma::uword i = 0; i < ref_depths.contigs.size(); ++i) {
    auto ref_ad = static_cast<s32>(ref_depths.obvs[i]);
    auto alt_ad = static_cast<s32>(alt_depths.obvs[i]);
    if (ref_ad + alt_ad == 0) {
      Logging::Debug("skipping variant at {}:{}-{} with zero total depth",
                     ref_depths.contigs[i],
                     ref_depths.starts[i],
                     ref_depths.ends[i]);
      continue;
    }
    auto baf = static_cast<f64>(alt_ad) / static_cast<f64>(ref_ad + alt_ad);
    f64 dh = 2 * std::abs(baf - 0.5);
    dh_obvs.contigs.push_back(ref_depths.contigs[i]);
    keep_starts.push_back(ref_depths.starts[i]);
    keep_ends.push_back(ref_depths.ends[i]);
    dh_vals.push_back(dh);
    dps.push_back(ref_ad + alt_ad);
  }
  dh_obvs.starts = arma::uvec(keep_starts);
  dh_obvs.ends = arma::uvec(keep_ends);
  dh_obvs.obvs = arma::vec(dh_vals);
  dh_obvs.dps = dps;
  dh_obvs.SetBAFSegObvsStatus(true);
  return dh_obvs;
}

/**
 * @brief read a VCF file, calculate Decrease in Heterozygosity (DH) statistics and return them in the form of a
 * Observations object
 * @param fname vcf file name
 * @param keep_all_variants keep homozygous variants
 * @param has_tumor_normal boolean to indicate if the VCF has a tumor and normal sample
 * @return Observations object
 */
Observations GetDhFromVcf(io::VcfReader& vcf_rdr,
                          bool keep_all_variants,
                          bool has_tumor_normal,
                          const BafFilterOptions& baf_filter_options,
                          const SampleMetadataOptions& sample_metadata_options) {
  Observations obvs;
  std::vector<arma::uword> starts;
  std::vector<arma::uword> ends;
  std::vector<f64> vals;
  auto update_fn = [&obvs, &starts, &ends, &vals](
                       const io::VcfRecordPtr& rec, const s32 ref_ad, const s32 alt_ad, const f64 baf) {
    obvs.contigs.push_back(rec->Chromosome());
    starts.push_back(rec->Position());
    ends.push_back(rec->Position() + 1);
    f64 dh = 2 * std::abs(baf - 0.5);
    vals.push_back(dh);
    obvs.dps.push_back(ref_ad + alt_ad);
  };

  if (has_tumor_normal) {
    if (!sample_metadata_options.normal_sample_name.has_value() ||
        !sample_metadata_options.tumor_sample_name.has_value()) {
      throw std::runtime_error("normal_sample_name and tumor_sample_name must be provided if in somatic T/N mode");
    }
    // VcfCheckResult checked by GetDepthsFromVcf; not needed here.
    (void)AccessBAFValuesFromTumorNormalVCF(vcf_rdr, update_fn, baf_filter_options, sample_metadata_options);
  } else {
    if (!sample_metadata_options.normal_sample_name.has_value()) {
      throw std::runtime_error("sample_name must be provided for single-sample VCF parsing");
    }
    AccessBAFValuesFromVCF(vcf_rdr, keep_all_variants, update_fn, baf_filter_options, sample_metadata_options);
  }
  obvs.starts = arma::uvec(starts);
  obvs.ends = arma::uvec(ends);
  obvs.obvs = arma::vec(vals);
  obvs.SetBAFSegObvsStatus(true);
  return obvs;
}

/**
 * @brief read a BED file containing reference and alternate allele depths for loci in the genome. No filtering is done
 * here. Expected BED format: Contig\tStart\tEnd\tBAF\tRef_AD\tAlt_AD (0-based half-open coordinates).
 * `#`-prefixed metadata/comments and a `#`-prefixed column header are supported and skipped.
 * @param fname BED file name
 * @return RefAltObservations object
 */
RefAltObservations GetDepthsFromFile(const fs::path& fname) {
  Observations ref_ads;
  Observations alt_ads;
  std::ifstream ifs(fname);
  if (!ifs.is_open()) {
    throw error::Error("Failed to open BAF BED file: {}", fname.string());
  }
  std::vector<arma::uword> starts;
  std::vector<arma::uword> ends;
  std::vector<f64> ref_vals;
  std::vector<f64> alt_vals;
  s32 row_num = 0;
  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream ss(line);
    std::string contig;
    std::string start_str;
    std::string end_str;
    std::string baf_str;
    std::string ref_str;
    std::string alt_str;
    if (!std::getline(ss, contig, '\t') || !std::getline(ss, start_str, '\t') || !std::getline(ss, end_str, '\t') ||
        !std::getline(ss, baf_str, '\t') || !std::getline(ss, ref_str, '\t') || !std::getline(ss, alt_str, '\t')) {
      throw error::Error("BAF BED file must have at least 6 columns (Contig, Start, End, BAF, Ref_AD, Alt_AD)");
    }
    // column 4 (BAF) is not read; it will be recalculated if needed
    ref_ads.contigs.push_back(contig);
    alt_ads.contigs.push_back(contig);
    try {
      starts.push_back(std::stoul(start_str));
      ends.push_back(std::stoul(end_str));
      ref_vals.push_back(std::stod(ref_str));
      alt_vals.push_back(std::stod(alt_str));
    } catch (const std::exception& e) {
      throw error::Error("BAF BED file parse error at row {}: {} (line: \"{}\")", row_num + 1, e.what(), line);
    }
    row_num += 1;
  }
  Logging::Debug("read {} variant records", row_num);
  ref_ads.starts = arma::uvec(starts);
  ref_ads.ends = arma::uvec(ends);
  ref_ads.obvs = arma::vec(ref_vals);
  ref_ads.SetBAFSegObvsStatus(true);
  alt_ads.starts = arma::uvec(starts);
  alt_ads.ends = arma::uvec(ends);
  alt_ads.obvs = arma::vec(alt_vals);
  alt_ads.SetBAFSegObvsStatus(true);
  return RefAltObservations{ref_ads, alt_ads, std::nullopt};
}

/**
 * @brief returns a new Observations object that represents the observations within the range [i, j) of the
 * current Observations object
 * @param i inclusive start of range
 * @param j exclusive end of range
 * @return a new Observations object containing only the values within the range [i,j)
 */
Observations Observations::FilterByRange(size_t i, size_t j) const {
  Observations ret;
  auto first = contigs.begin();
  auto last = contigs.begin();
  std::advance(first, i);
  std::advance(last, j);
  ret.contigs = std::vector(first, last);
  ret.starts = starts(arma::span(i, j - 1));
  ret.ends = ends(arma::span(i, j - 1));
  ret.obvs = obvs(arma::span(i, j - 1));
  return ret;
}

/**
 * @brief returns a new Observations object that represents a subset of the observations from the given
 * indices in current Observations object
 * @param idxs - list of indexes to extract into a new object
 * @return a new Observations object containing only the values from the original object represented by idxs
 */
Observations Observations::FilterByIdxs(const std::vector<size_t>& idxs) const {
  Observations ret;
  ret.contigs.resize(idxs.size());
  ret.starts.resize(idxs.size());
  ret.ends.resize(idxs.size());
  ret.obvs.resize(idxs.size());
  size_t i = 0;
  for (const auto& j : idxs) {
    ret.contigs[i] = contigs[j];
    ret.starts[i] = starts[j];
    ret.ends[i] = ends[j];
    ret.obvs[i] = obvs[j];
    i += 1;
  }
  return ret;
}

/**
 * @brief Write the observations to a BED file with `##`-prefixed metadata and `#`-prefixed column header.
 * @param ofs output file stream
 * @param value_column_name name of the value column
 * @param value_as_int whether to write the value as an integer
 * @param command_line_info metadata to write at the top of the file
 */
void Observations::Write(std::ofstream& ofs,
                         const std::string& value_column_name,
                         const bool value_as_int,
                         const io::CommandLineInfo& command_line_info) const {
  auto writer = csv::make_tsv_writer(ofs);
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info);
  }
  vec<std::string> column_names{kColumnContigAsComment, kColumnStart, kColumnEnd, value_column_name};
  writer << column_names;
  for (size_t i = 0; i < contigs.size(); ++i) {
    std::string value_string = value_as_int ? FloatAsIntString(obvs[i]) : FloatAsString(obvs[i]);
    writer << std::make_tuple(contigs[i], starts[i], ends[i], value_string);
  }
}

bool Observations::IsSorted() const {
  std::set<std::string> seen_contigs;
  for (size_t i = 1; i < contigs.size(); ++i) {
    if (contigs[i] != contigs[i - 1]) {
      seen_contigs.insert(contigs[i - 1]);
    }
    if (seen_contigs.find(contigs[i]) != seen_contigs.end()) {
      Logging::Warn("contigs not ordered");
      return false;
    }
    if (starts[i] <= starts[i - 1]) {
      Logging::Warn("starts not ordered");
      return false;
    }
    if (ends[i] <= ends[i - 1]) {
      Logging::Warn("ends not ordered");
      return false;
    }
  }
  return true;
}

void Observations::PopulateFieldsFromRegions() {
  contigs.resize(regions.size());
  starts.resize(regions.size());
  ends.resize(regions.size());
  for (size_t i = 0; i < regions.size(); ++i) {
    auto [contig, start, end] = ParseRegionString(regions[i]);
    contigs[i] = contig;
    starts[i] = start;
    ends[i] = end;
  }
}

/**
 * @brief Extract tumor and normal sample information from VCF header
 * @param header Shared pointer to VcfHeader
 * @return VcfHeaderInfo struct containing sample indexes for `tumor_sample` and `normal_sample` in the VCF file.
 */
VcfHeaderInfo GetTumorNormalSampleInfoFromHeader(const std::shared_ptr<io::VcfHeader>& header,
                                                 const std::string& normal_sample_name,
                                                 const std::string& tumor_sample_name) {
  VcfHeaderInfo info{};
  // extract sample indexes for `tumor_sample` and `normal_sample` if they are found in the header
  const std::map<std::string, s32>& sample_indexes{header->GetSampleIndexes()};
  if (sample_indexes.find(normal_sample_name) == sample_indexes.end()) {
    Logging::Error("Invalid VCF sample index for normal_sample ({})", normal_sample_name);
    throw std::runtime_error("Invalid VCF sample index");
  }
  if (sample_indexes.find(tumor_sample_name) == sample_indexes.end()) {
    Logging::Error("Invalid VCF sample index for tumor_sample ({})", tumor_sample_name);
    throw std::runtime_error("Invalid VCF sample index");
  }
  info.normal_index = sample_indexes.at(normal_sample_name);
  info.tumor_index = sample_indexes.at(tumor_sample_name);
  if (info.tumor_index == info.normal_index) {
    Logging::Error("tumor_sample ({}) and normal_sample ({}) have the same VCF sample index: {}",
                   tumor_sample_name,
                   normal_sample_name,
                   info.normal_index);
    throw std::runtime_error("Invalid VCF sample index");
  }
  info.has_tumor_normal = true;
  Logging::Info("normal_sample: {}", normal_sample_name);
  Logging::Info("tumor_sample:  {}", tumor_sample_name);
  return info;
}

/**
 * @brief Extract single-sample information from a VCF header, verifying that the requested sample name exists
 * @param header Shared pointer to VcfHeader
 * @param normal_sample_name Name of the sample to look up
 * @return VcfHeaderInfo struct with normal_index set to the column index of the requested sample
 */
VcfHeaderInfo GetSingleSampleInfoFromHeader(const std::shared_ptr<io::VcfHeader>& header,
                                            const std::string& normal_sample_name) {
  VcfHeaderInfo info{};
  const std::map<std::string, s32>& sample_indexes{header->GetSampleIndexes()};
  if (sample_indexes.find(normal_sample_name) == sample_indexes.end()) {
    Logging::Error("Invalid VCF sample index for normal_sample ({})", normal_sample_name);
    throw std::runtime_error("Invalid VCF sample index");
  }
  info.normal_index = sample_indexes.at(normal_sample_name);
  info.normal_sample = normal_sample_name;
  info.has_tumor_normal = false;
  Logging::Info("normal_sample: {}", normal_sample_name);
  return info;
}

/**
 * @brief Write the observations to a BigWig file
 * @param filename path to output BigWig file
 * @param fai_filename path to FASTA index file containing chromosome lengths
 */
void Observations::WriteBigWig(const fs::path& filename, const fs::path& fai_filename) const {
  ::xoos::cnc::WriteObservationsToBigWig(*this, filename, fai_filename);
}

f64 ComputeExtremeBafProportion(const Observations& ref_depths, const Observations& alt_depths) {
  if (ref_depths.contigs.size() != alt_depths.contigs.size()) {
    throw error::Error(
        "Reference allele depths and alternate allele depths must have the same number of observations "
        "(reference: {}, alternate: {})",
        ref_depths.contigs.size(),
        alt_depths.contigs.size());
  }
  if (ref_depths.contigs.size() != ref_depths.obvs.n_elem || alt_depths.contigs.size() != alt_depths.obvs.n_elem) {
    throw error::Error(
        "Allele depth observations are internally inconsistent: contig count does not match observation count "
        "(reference contigs: {}, reference observations: {}, alternate contigs: {}, alternate observations: {})",
        ref_depths.contigs.size(),
        ref_depths.obvs.n_elem,
        alt_depths.contigs.size(),
        alt_depths.obvs.n_elem);
  }
  if (ref_depths.contigs.empty()) {
    return 0.0;
  }
  size_t extreme_count = 0;
  size_t valid_count = 0;
  for (arma::uword i = 0; i < ref_depths.contigs.size(); ++i) {
    const auto ref_ad = static_cast<f64>(ref_depths.obvs[i]);
    const auto alt_ad = static_cast<f64>(alt_depths.obvs[i]);
    const f64 total = ref_ad + alt_ad;
    if (total <= 0.0) {
      continue;
    }
    ++valid_count;
    const f64 baf = alt_ad / total;
    if (baf > kExtremeBafHighThreshold || baf < kExtremeBafLowThreshold) {
      ++extreme_count;
    }
  }
  if (valid_count == 0) {
    return 0.0;
  }
  return static_cast<f64>(extreme_count) / static_cast<f64>(valid_count);
}

bool CheckExtremeBafProportionAndWarn(const Observations& ref_depths,
                                      const Observations& alt_depths,
                                      const f64 purity) {
  // At high purity, extreme BAF values are biologically expected due to LOH
  // so the extreme proportion check is not meaningful — skip the warning.
  if (purity >= kExtremeBafPurityThreshold) {
    return false;
  }
  const f64 extreme_baf_proportion = ComputeExtremeBafProportion(ref_depths, alt_depths);
  // A low extreme BAF proportion is normal; only warn when it exceeds the threshold,
  // which may indicate a tumor-normal sample mismatch.
  if (extreme_baf_proportion <= kExtremeBafProportionThreshold) {
    return false;
  }
  Logging::Warn(
      "Sample has a high proportion ({:.1f}%) of heterozygous SNPs with extreme tumor BAF (BAF > {:.2f} or BAF < "
      "{:.2f}) at predicted purity ({:.2f}). The tumor and normal samples may not be from the same individual.",
      extreme_baf_proportion * 100.0,
      kExtremeBafHighThreshold,
      kExtremeBafLowThreshold,
      purity);
  return true;
}
}  // namespace xoos::cnc
