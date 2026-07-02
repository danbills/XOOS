#include "calculate-coverage/calculate-coverage.h"

#include <fstream>
#include <limits>

#include <htslib/sam.h>

#include <csv.hpp>

#include <taskflow/algorithm/for_each.hpp>

#include <xoos/error/error.h>
#include <xoos/io/alignment-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/util/file-functions.h>

#include "baits.h"
#include "calculate-coverage/calculate-coverage-options.h"
#include "calculate-coverage/calculate-coverage-taskflow-graph.h"
#include "copy-number-caller/copy-number-caller-options.h"
#include "io/column-names.h"
#include "io/copy-number-caller-default-filenames.h"
#include "utility/utility-functions.h"

namespace xoos::cnc {

/// Returns true if the record has a YF:i:1 BAM auxiliary tag, indicating a rescued
/// unprojectable read from pangenome alignment.
static bool HasRescuedReadTag(const bam1_t* record) {
  const auto* tag_ptr = bam_aux_get(record, "YF");
  return tag_ptr != nullptr && bam_aux2i(tag_ptr) == 1;
}

static s32 ReadBam(void* data, bam1_t* b) {
  auto* plp = static_cast<Plpconf*>(data);
  s32 ret;
  u32 dup_flag = 1024;
  bool keep_reading = true;
  while (keep_reading) {
    ret = sam_itr_next(plp->fp, plp->itr, b);
    if (ret < 0) {
      // actual read error (not EOF) — record it so the caller can throw
      if (ret < -1) {
        plp->had_read_error = true;
      }
      break;
    }
    if (b->core.flag & dup_flag) {
      plp->marked_dup = true;
    }
    bool excluded = (b->core.flag & plp->exclude_flags) != 0;
    // Rescue secondary reads tagged with YF:i:1 (unprojectable pangenome reads).
    // Re-apply remaining exclude flags so rescued reads are still filtered for
    // duplicates, unmapped, QC-fail, etc. — same as primary reads.
    if (excluded && (plp->exclude_flags & BAM_FSECONDARY) != 0 && (b->core.flag & BAM_FSECONDARY) != 0 &&
        HasRescuedReadTag(b)) {
      const u32 exclude_without_secondary = plp->exclude_flags & ~static_cast<u32>(BAM_FSECONDARY);
      excluded = (b->core.flag & exclude_without_secondary) != 0;
      if (!excluded) {
        plp->has_rescued_secondaries = true;
      }
    }
    keep_reading = excluded;
  }
  return ret;
}

static s32 IncreaseReadCount(void* data, const bam1_t* b, bam_pileup_cd* cd) {
  auto* plp = static_cast<Plpconf*>(data);
  plp->total_reads++;
  plp->sum_mapping_quality += b->core.qual;
  return 0;
}

static void WriteMapqs(const CoverageRecords& cov,
                       const io::CommandLineInfo& command_line_info,
                       const fs::path& out_file) {
  std::ofstream ofs(out_file);
  if (!ofs.is_open()) {
    Logging::Error("Failed to open mapping quality output file: {}", out_file.string());
    throw std::runtime_error("Failed to open mapping quality output file");
  }
  auto writer = csv::make_tsv_writer(ofs);
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info);
  }
  vec<std::string> column_names{kColumnContigAsComment, kColumnStart, kColumnEnd, kColumnMeanMapq};
  writer << column_names;
  for (size_t i = 0; i < cov.region.size(); ++i) {
    const auto& [contig, start, end] = ParseRegionString(cov.region[i]);
    writer << std::make_tuple(contig, start, end, cov.mean_mapping_quality[i]);
  }
}

void CalculateCoverageMain(const CopyNumberCallerOptions& options) {
  Logging::Info("Reading baits file");
  std::ifstream bait_istream(options.augmented_baits_fname.value());
  BaitRecords baits(bait_istream);

  const auto coverage_out = options.output_dir / kDefaultCoverageOutput;
  const auto mapping_qualities_out = options.output_dir / kDefaultMapQOutput;
  file::CheckFilePermissionsAndOutputPathExistence({}, {coverage_out, mapping_qualities_out});

  tf::Taskflow taskflow;
  CalculateCoverageTaskflowGraph calculate_coverage_tf_graph(options.normal_bam_fname,
                                                             baits,
                                                             options.calculate_coverage_options.exclude_flags,
                                                             options.calculate_coverage_options.ignore_DN);
  taskflow.composed_of(calculate_coverage_tf_graph);
  tf::Executor executor(options.threads);
  executor.run(taskflow).get();
  Logging::Info("Writing coverages to {}", coverage_out.string());
  CoverageRecords cov = calculate_coverage_tf_graph.GetResult();
  if (!cov.has_duplicates) {
    Logging::Info("The input BAM file does not have any reads marked as duplicate");
  }
  std::ofstream ofs(coverage_out);
  if (!ofs.is_open()) {
    throw error::Error("Failed to open coverage output file: {}", coverage_out.string());
  }
  WriteMapqs(cov, options.command_line_info, mapping_qualities_out);
  cov.Write(ofs, options.command_line_info);
}

static s32 CountSkippedBases(const bam_pileup1_t* plp_column, const s32 n) {
  s32 n_skipped = 0;
  for (s32 i = 0; i < n; ++i) {
    if (plp_column[i].is_del == 1 || plp_column[i].is_refskip == 1) {
      n_skipped++;
    }
  }
  return n_skipped;
}

static s32 CalculateSingleRegionCoverage(Plpconf& plp, const s64 start, const s64 end, const bool ignore_dn) {
  io::BamPlpPtr plp_itr(bam_plp_init(ReadBam, static_cast<void*>(&plp)));
  bam_plp_set_maxcnt(plp_itr.get(),
                     INT_MAX);  // Note: default maxcnt is 8000, but set to INT_MAX to mimic samtools bedcov behavior
  bam_plp_constructor(plp_itr.get(), IncreaseReadCount);

  s32 tid = -1;
  s32 pos = -1;
  s32 n = 0;
  s32 total_bases = 0;
  const bam_pileup1_t* plp_column;

  while ((plp_column = bam_plp_auto(plp_itr.get(), &tid, &pos, &n))) {
    if (pos >= start && pos < end) {
      const s32 n_skipped = ignore_dn ? CountSkippedBases(plp_column, n) : 0;
      total_bases += (n - n_skipped);
    }
  }
  return total_bases;
}

void CalculateCoverageRegion(const io::AlignmentReader& reader,
                             CoverageRecords& res,
                             const u32 exclude_flags,
                             const bool ignore_dn,
                             const s32 start_idx,
                             const s32 step) {
  Plpconf plp;
  plp.fp = reader.bam.get();
  plp.marked_dup = false;

  const auto total_regions = static_cast<s32>(res.region.size());
  const s32 end_idx = std::min(start_idx + step, total_regions);

  for (s32 i = start_idx; i < end_idx; ++i) {
    const auto idx = static_cast<arma::uword>(i);
    // baits file is 1-based start, and gets converted to 0-based half-open with ParseRegionString
    auto [chr, start_raw, end_raw] = ParseRegionString(res.region[idx]);
    const auto start = static_cast<s64>(start_raw);
    const auto end = static_cast<s64>(end_raw);

    if (end <= start) {
      throw std::runtime_error("Error: record end < start in regions file");
    }
    const auto tid = io::SamHdrName2Tid(reader.hdr.get(), chr);
    if (tid < 0) {
      Logging::Warn("Contig {} not found in bam!", chr);
    }

    const auto itr = io::SamItrQueryI(reader.idx.get(), tid, start, end);
    plp.itr = itr.get();
    plp.exclude_flags = exclude_flags;
    plp.total_reads = 0;
    plp.sum_mapping_quality = 0;
    plp.has_rescued_secondaries = false;

    const s32 total_bases = CalculateSingleRegionCoverage(plp, start, end, ignore_dn);

    if (plp.had_read_error) {
      throw error::Error("Error reading BAM file '{}': file may be truncated or corrupted", reader.bam->fn);
    }
    if (plp.marked_dup) {
      res.has_duplicates = true;
    }
    if (plp.has_rescued_secondaries) {
      res.has_rescued_secondaries = true;
    }
    res.total_coverage[idx] = total_bases;
    res.count[idx] = plp.total_reads;
    res.mean_mapping_quality[idx] = plp.total_reads > 0
                                        ? static_cast<f64>(plp.sum_mapping_quality) / static_cast<f64>(plp.total_reads)
                                        : std::numeric_limits<f64>::quiet_NaN();
  }
}

}  // namespace xoos::cnc
