#include "io/write-metrics.h"

#include <fstream>

#include <fmt/core.h>

#include <csv.hpp>

#include <xoos/log/logging.h>

#include "likelihood/likelihood-flags.h"
#include "seg-to-vcf/seg-to-vcf.h"

namespace xoos::cnc {

void WriteMetricsFile(const fs::path& metrics_fname,
                      const BaitRecords& baits,
                      size_t n_snps,
                      Sex sex,
                      const std::vector<GenomicSegment>& segments,
                      const std::optional<f64>& purity,
                      const std::optional<f64>& ploidy,
                      CopyNumberCallerModes mode,
                      const io::CommandLineInfo& command_line_info,
                      const std::optional<CoverageCheckResult>& coverage_check,
                      const std::optional<VcfCheckResult>& vcf_check) {
  std::ofstream metrics_ofs(metrics_fname);
  if (!metrics_ofs.is_open()) {
    Logging::Error("Failed to open metrics file: {}", metrics_fname.string());
    throw std::runtime_error("Failed to open metrics file: " + metrics_fname.string());
  }
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(metrics_ofs, command_line_info);
  }
  std::size_t num_gain_segments = 0;
  std::size_t num_loss_segments = 0;
  std::size_t num_cnloh_segments = 0;
  std::size_t num_gainloh_segments = 0;
  std::size_t num_passing_gain_segments = 0;
  std::size_t num_passing_loss_segments = 0;
  std::size_t num_passing_cnloh_segments = 0;
  std::size_t num_passing_gainloh_segments = 0;
  segmentation::CopyNumberStatus status = segmentation::CopyNumberStatus::kUnknown;
  for (const auto& seg : segments) {
    if (mode == CopyNumberCallerModes::kGermlineNormalWGS) {
      status = segmentation::GetCopyNumberStatus(seg);
    } else if (mode == CopyNumberCallerModes::kSomaticTumorNormalWGS) {
      status = segmentation::GetSomaticCopyNumberStatus(seg);
    }
    bool is_pass = segmentation::ConstructVcfFilterFieldFromGenomicSegment(seg) == kLikelihoodFlagPass;
    if (status == segmentation::CopyNumberStatus::kGain) {
      ++num_gain_segments;
      num_passing_gain_segments += is_pass ? 1 : 0;
    } else if (status == segmentation::CopyNumberStatus::kLoss) {
      ++num_loss_segments;
      num_passing_loss_segments += is_pass ? 1 : 0;
    } else if (status == segmentation::CopyNumberStatus::kCnLoh) {
      ++num_cnloh_segments;
      num_passing_cnloh_segments += is_pass ? 1 : 0;
    } else if (status == segmentation::CopyNumberStatus::kGainLoh) {
      ++num_gainloh_segments;
      num_passing_gainloh_segments += is_pass ? 1 : 0;
    }
  }

  if (mode == CopyNumberCallerModes::kGermlineNormalWGS) {
    metrics_ofs << "Metric\tValue\n"
                << "Number of targets\t" << baits.GetRegions().size() << '\n'
                << "Sex\t" << SexToStr(sex) << '\n'
                << "Number of segments\t" << segments.size() << '\n'
                << "Number of gain segments\t" << num_gain_segments << '\n'
                << "Number of loss segments\t" << num_loss_segments << '\n'
                << "Number of passing gain segments\t" << num_passing_gain_segments << '\n'
                << "Number of passing loss segments\t" << num_passing_loss_segments << '\n';
  } else if (mode == CopyNumberCallerModes::kSomaticTumorNormalWGS) {
    const std::string purity_str = purity.has_value() ? fmt::format("{:.2f}", purity.value()) : "NA";
    const std::string ploidy_str = ploidy.has_value() ? fmt::format("{:.2f}", ploidy.value()) : "NA";
    metrics_ofs << "Metric\tValue\n"
                << "Number of targets\t" << baits.GetRegions().size() << '\n'
                << "Number of SNPs\t" << n_snps << '\n'
                << "Sex\t" << SexToStr(sex) << '\n'
                << "Tumor Purity\t" << purity_str << '\n'
                << "Tumor Ploidy\t" << ploidy_str << '\n'
                << "Number of segments\t" << segments.size() << '\n'
                << "Number of somatic gain (non-LOH) segments\t" << num_gain_segments << '\n'
                << "Number of somatic loss segments\t" << num_loss_segments << '\n'
                << "Number of somatic copy-neutral LOH segments\t" << num_cnloh_segments << '\n'
                << "Number of somatic gain LOH segments\t" << num_gainloh_segments << '\n'
                << "Number of passing somatic gain (non-LOH) segments\t" << num_passing_gain_segments << '\n'
                << "Number of passing somatic loss segments\t" << num_passing_loss_segments << '\n'
                << "Number of passing somatic copy-neutral LOH segments\t" << num_passing_cnloh_segments << '\n'
                << "Number of passing somatic gain LOH segments\t" << num_passing_gainloh_segments << '\n';
  }

  if (coverage_check.has_value()) {
    metrics_ofs << fmt::format("Median coverage\t{:.2f}\n", coverage_check->median_coverage)
                << fmt::format("Percentage of low coverage windows\t{:.2f}\n",
                               coverage_check->pct_windows_below_threshold);
  }

  if (vcf_check.has_value()) {
    const auto& vc = *vcf_check;
    const std::string het_snp_fraction_str =
        vc.het_snp_fraction.has_value() ? fmt::format("{:.4f}", *vc.het_snp_fraction) : "NA";
    metrics_ofs << "Number of heterozygous germline SNPs\t" << vc.num_het_snps << '\n'
                << "Number of somatic variants\t" << vc.num_somatic_variants << '\n'
                << "Het SNP fraction\t" << het_snp_fraction_str << '\n';
  }
}

}  // namespace xoos::cnc
