#pragma once
#include <optional>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>

#include "copy-number-caller/common/coverage-check.h"
#include "copy-number-caller/common/vcf-check.h"
#include "copy-number-caller/germline-normal-wgs/germline-normal-wgs.h"
#include "copy-number-caller/somatic-tumor-normal-wgs/somatic-tumor-normal-wgs.h"

namespace xoos::cnc {

/**
 * @brief Write all output files for a germline workflow early exit.
 *
 * Writes header-only versions of all expected output files so downstream
 * consumers don't fail on missing files. The coverage file is already written
 * before this is called.
 */
void WriteGermlineEmptyOutputs(const GermlineNormalWGSOutputPaths& paths,
                               const CoverageCheckResult& coverage_check,
                               const BaitRecords& baits,
                               const CopyNumberCallerOptions& options);

/**
 * @brief Write all output files for a somatic workflow early exit.
 *
 * Handles both low-coverage and VCF data-quality graceful exits.
 * When coverage_check is provided, coverage metrics are included in the metrics file.
 * When vcf_check is provided, VCF-related metrics are included in the metrics file.
 */
void WriteSomaticEmptyOutputs(const SomaticTumorNormalWGSOutputPaths& paths,
                              const std::optional<CoverageCheckResult>& coverage_check,
                              const BaitRecords& baits,
                              const CopyNumberCallerOptions& options,
                              const std::optional<VcfCheckResult>& vcf_check = std::nullopt);

}  // namespace xoos::cnc
