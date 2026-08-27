#pragma once

/**
 * ============================================================================
 * Copy Number Caller & HMM Segmentation Core Types
 * File: copy_number_caller/src/cuda/cnv_types.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Design Goals:
 *   - 16-byte aligned genomic bin records (Log2R, GC content, BAF, weights).
 *   - Discrete Copy Number HMM states (0=HOM_DEL, 1=HET_DEL, 2=DIPLOID, 3=GAIN, 4=AMP, 5=HIGH_AMP).
 *   - 2D grid evaluation structures for tumor purity (p) and ploidy (psi).
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <cuda_runtime.h>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace xoos::cnv_caller::cuda {

constexpr uint32_t kMaxHMMStates = 8;        // Copy number states 0, 1, 2, 3, 4, 5, 6, 7
constexpr uint32_t kPurityGridSteps = 100;    // Purity scan: 0.01 to 1.00 (step 0.01)
constexpr uint32_t kPloidyGridSteps = 100;    // Ploidy scan: 1.00 to 6.00 (step 0.05)

/**
 * @brief Discrete Copy Number States.
 */
enum class CopyNumberState : uint8_t {
    HOMOZYGOUS_DELETION = 0, // CN = 0
    HETEROZYGOUS_DELETION = 1, // CN = 1
    NEUTRAL_DIPLOID = 2,     // CN = 2
    SINGLE_COPY_GAIN = 3,    // CN = 3
    AMPLIFICATION = 4,       // CN = 4
    HIGH_AMPLIFICATION = 5   // CN >= 5
};

/**
 * @brief Genomic Bin Record for CNV analysis.
 */
struct __align__(16) GenomicBinRecord {
    uint32_t chr_id;
    uint32_t start_pos;
    uint32_t end_pos;
    float gc_content;       // [0.0, 1.0]
    float raw_depth;        // Observed read count or coverage
    float normalized_log2r; // GC-corrected and PoN normalized Log2 ratio
    float baf;              // B-allele frequency [0.0, 0.5] (mirror folded)
    float weight;           // Inverse variance weight
};

/**
 * @brief Segment Output Record.
 */
struct __align__(16) CnvSegment {
    uint32_t chr_id;
    uint32_t start_pos;
    uint32_t end_pos;
    uint32_t num_bins;
    float mean_log2r;
    float mean_baf;
    uint8_t copy_number;    // Inferred integer copy number
    float confidence;       // Posterior probability or Phred quality
};

/**
 * @brief Tumor Purity & Ploidy Global Fit.
 */
struct __align__(16) PurityPloidyFit {
    float best_purity;      // Estimated tumor purity (0.0 to 1.0)
    float best_ploidy;      // Estimated tumor ploidy (e.g. 2.0 = diploid)
    double max_log_likelihood;
    float subclonal_fraction;
};

/**
 * @brief CNV Execution Statistics.
 */
struct CnvExecutionStats {
    uint64_t total_bins_processed = 0;
    uint64_t total_segments_discovered = 0;
    uint64_t total_grid_evaluations = 0;
    double gc_correction_time_ms = 0.0;
    double hmm_viterbi_time_ms = 0.0;
    double purity_ploidy_time_ms = 0.0;
    double total_time_ms = 0.0;
    double throughput_bins_per_sec = 0.0;
};

/**
 * @brief Host utility: Collapse contiguous state sequence into segment descriptors.
 */
inline void collapse_state_path_to_segments(
    const GenomicBinRecord* h_bins,
    const uint8_t* h_state_path,
    uint64_t num_bins,
    std::vector<CnvSegment>& out_segments
) {
    out_segments.clear();
    if (num_bins == 0) return;

    size_t seg_start = 0;
    uint8_t current_cn = h_state_path[0];
    uint32_t current_chr = h_bins[0].chr_id;

    for (size_t i = 1; i <= num_bins; ++i) {
        bool boundary = (i == num_bins) ||
                        (h_state_path[i] != current_cn) ||
                        (h_bins[i].chr_id != current_chr);

        if (boundary) {
            size_t seg_len = i - seg_start;
            double sum_log2r = 0.0;
            double sum_baf = 0.0;

            for (size_t k = seg_start; k < i; ++k) {
                sum_log2r += h_bins[k].normalized_log2r;
                sum_baf += h_bins[k].baf;
            }

            CnvSegment seg;
            seg.chr_id = current_chr;
            seg.start_pos = h_bins[seg_start].start_pos;
            seg.end_pos = h_bins[i - 1].end_pos;
            seg.num_bins = static_cast<uint32_t>(seg_len);
            seg.mean_log2r = static_cast<float>(sum_log2r / seg_len);
            seg.mean_baf = static_cast<float>(sum_baf / seg_len);
            seg.copy_number = current_cn;
            seg.confidence = 99.0f;

            out_segments.push_back(seg);

            if (i < num_bins) {
                seg_start = i;
                current_cn = h_state_path[i];
                current_chr = h_bins[i].chr_id;
            }
        }
    }
}

} // namespace xoos::cnv_caller::cuda
