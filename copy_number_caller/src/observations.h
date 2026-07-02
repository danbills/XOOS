#pragma once
#include <armadillo>
#include <memory>
#include <string>
#include <vector>

#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "copy-number-caller/common/vcf-check.h"
#include "misc/sample-metadata-options.h"
#include "misc/vcf-parsing-options.h"

namespace xoos::cnc {

struct Observations {
 public:
  Observations() : _is_baf_seg_obs(false) {
  }

  Observations(std::vector<std::string>&& c, const arma::uvec&& s, const arma::uvec&& e, const arma::vec&& o)
      : contigs(c), regions({}), starts(s), ends(e), obvs(o), dps({}), _is_baf_seg_obs(false) {
  }

  explicit Observations(size_t size) {
    contigs.resize(size);
    regions.resize(size);
    starts.resize(size);
    ends.resize(size);
    obvs.resize(size);
    dps.resize(size);
    _is_baf_seg_obs = false;
  }

  Observations(const Observations& rhs) = default;
  Observations(Observations&& rhs) noexcept = default;
  Observations& operator=(const Observations& rhs) = default;
  Observations& operator=(Observations&& rhs) noexcept = default;
  ~Observations() = default;

  void SetBAFSegObvsStatus(bool status) {
    _is_baf_seg_obs = status;
  }

  void PopulateFieldsFromRegions();

  bool IsBAFSegObvs() const {
    return _is_baf_seg_obs;
  }

  bool IsSorted() const;

  std::vector<std::string> contigs;
  std::vector<std::string> regions;
  arma::uvec starts;
  arma::uvec ends;
  arma::vec obvs;
  std::vector<s32> dps;
  Observations FilterByRange(size_t i, size_t j) const;
  Observations FilterByIdxs(const std::vector<size_t>& idxs) const;

  void Write(std::ofstream& ofs,
             const std::string& value_column_name,
             bool value_as_int,
             const io::CommandLineInfo& command_line_info) const;
  void WriteBigWig(const fs::path& filename, const fs::path& fai_filename) const;

 private:
  bool _is_baf_seg_obs;
};

struct RefAltObservations {
  Observations ref_obvs;
  Observations alt_obvs;

  // optional vector of somatic VAFs, only populated if the VCF is from a tumor-normal pair and contains somatic
  // variants.
  std::optional<std::vector<f64>> somatic_vafs;

  // VCF data-quality check result, populated for tumor-normal VCF parsing.
  std::optional<VcfCheckResult> vcf_check;
};

struct VcfHeaderInfo {
  std::string normal_sample;
  std::string tumour_sample;
  bool has_tumor_normal{false};
  s32 normal_index{-1};
  s32 tumor_index{-1};
};

Observations ReadObservations(std::istream& in_stream);
Observations GetDhFromDepths(const Observations& ref_depths, const Observations& alt_depths);
Observations GetDhFromVcf(io::VcfReader& vcf_rdr,
                          bool keep_all_variants,
                          bool has_tumor_normal,
                          const BafFilterOptions& baf_filter_options,
                          const SampleMetadataOptions& sample_metadata_options);
RefAltObservations GetDepthsFromVcf(io::VcfReader& vcf_rdr,
                                    bool keep_all_variants,
                                    bool has_tumor_normal,
                                    const BafFilterOptions& baf_filter_options,
                                    const SampleMetadataOptions& sample_metadata_options);

VcfHeaderInfo GetTumorNormalSampleInfoFromHeader(const std::shared_ptr<io::VcfHeader>& header,
                                                 const std::string& normal_sample_name,
                                                 const std::string& tumor_sample_name);
VcfHeaderInfo GetSingleSampleInfoFromHeader(const std::shared_ptr<io::VcfHeader>& header,
                                            const std::string& normal_sample_name);
RefAltObservations GetDepthsFromFile(const fs::path& fname);

constexpr f64 kExtremeBafLowThreshold = 0.03;
constexpr f64 kExtremeBafHighThreshold = 0.97;
constexpr f64 kExtremeBafProportionThreshold = 0.10;
constexpr f64 kExtremeBafPurityThreshold = 0.80;

/*
 * @brief Compute the proportion of heterozygous SNPs (based on normal) with extreme tumor BAF (> 0.97 or < 0.03).
 * @param ref_depths Reference allele depths from the tumor
 * @param alt_depths Alternate allele depths from the tumor
 * @return Proportion of SNPs with extreme BAF, or 0.0 if there are no valid depth loci
 * @throws error::Error if ref_depths and alt_depths have mismatched sizes
 */
f64 ComputeExtremeBafProportion(const Observations& ref_depths, const Observations& alt_depths);

/*
 * @brief Check if the extreme BAF proportion is high and log a warning if so.
 * @param ref_depths Reference allele depths from the tumor
 * @param alt_depths Alternate allele depths from the tumor
 * @param purity Estimated tumor purity
 * @return True if the extreme BAF proportion is high, false otherwise
 */
bool CheckExtremeBafProportionAndWarn(const Observations& ref_depths, const Observations& alt_depths, f64 purity);
}  // namespace xoos::cnc
