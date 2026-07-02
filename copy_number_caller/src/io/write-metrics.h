#pragma once
#include <optional>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>

#include "baits.h"
#include "copy-number-caller/common/coverage-check.h"
#include "copy-number-caller/common/vcf-check.h"
#include "copy-number-caller/copy-number-caller-modes.h"
#include "segmentation/genomic-segments.h"
#include "sex.h"

namespace xoos::cnc {
using segmentation::GenomicSegment;
/**
 * @brief Writes summary metrics about the copy number run
 *
 * @param metrics_fname fs::path to the output metrics file.
 * @param baits BaitRecords object containing target regions.
 * @param n_snps number of SNPs used for analysis
 * @param sex Sample sex
 * @param segments Vector of GenomicSegment objects to analyze.
 * @param purity Sample purity. Only used for somatic modes.
 * @param ploidy Sample ploidy. Only used for somatic modes.
 * @param mode The mode of the CopyNumberCaller
 * @param command_line_info Metadata to write at the top of the file.
 * @param coverage_check Optional coverage check result. When provided, appends
 *        "Median coverage" and "Percentage of low coverage windows" rows.
 * @param vcf_check Optional VCF check result. When provided, appends
 *        "Number of heterozygous germline SNPs", "Number of somatic variants",
 *        and "Het SNP fraction" rows. Only applicable for tumor-normal mode.
 */
void WriteMetricsFile(const fs::path& metrics_fname,
                      const BaitRecords& baits,
                      size_t n_snps,
                      Sex sex,
                      const std::vector<GenomicSegment>& segments,
                      const std::optional<f64>& purity,
                      const std::optional<f64>& ploidy,
                      CopyNumberCallerModes mode,
                      const io::CommandLineInfo& command_line_info,
                      const std::optional<CoverageCheckResult>& coverage_check = std::nullopt,
                      const std::optional<VcfCheckResult>& vcf_check = std::nullopt);
}  // namespace xoos::cnc
