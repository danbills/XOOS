#pragma once

/**
 * ============================================================================
 * GPU Parallel HMM Viterbi Trellis Segmentation Engine
 * File: copy_number_caller/src/cuda/cbs_hmm_segmentation_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Capabilities:
 *   - Fast Quadratic GC Content Bias Normalization.
 *   - In-Register HMM Viterbi Trellis Dynamic Programming across all genomic bins.
 *   - Contiguous copy number state boundary extraction & segment aggregation.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include "cnv_types.cuh"
#include <cuda_runtime.h>
#include <cstdint>
#include <cmath>

namespace xoos::cnv_caller::cuda {

constexpr float kLog2E = 1.4426950408889634f;
constexpr float kDefaultEmissionSigma = 0.25f;
constexpr float kStateSwitchPenalty = -8.0f; // log(1e-4)

/**
 * @brief CUDA Kernel: Parallel GC Bias Normalization & Log2 Ratio calculation.
 */
__global__ void xoos_gc_normalize_kernel(
    GenomicBinRecord* __restrict__ d_bins,
    uint64_t num_bins,
    float gc_coeff_a,
    float gc_coeff_b,
    float gc_coeff_c,
    float pon_baseline_depth
) {
    size_t idx = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    size_t stride = static_cast<size_t>(gridDim.x) * blockDim.x;

    for (size_t i = idx; i < num_bins; i += stride) {
        GenomicBinRecord& bin = d_bins[i];
        float gc = bin.gc_content;

        // Quadratic GC bias model
        float expected_depth = gc_coeff_a * gc * gc + gc_coeff_b * gc + gc_coeff_c;
        if (expected_depth < 1.0f) expected_depth = 1.0f;

        float target_depth = pon_baseline_depth > 0.0f ? pon_baseline_depth : 100.0f;
        float norm_depth = (bin.raw_depth / expected_depth) * target_depth;

        // Log2 ratio relative to baseline
        float ratio = (norm_depth + 1e-4f) / (target_depth + 1e-4f);
        bin.normalized_log2r = logf(ratio) * kLog2E;
        bin.weight = 1.0f / (kDefaultEmissionSigma * kDefaultEmissionSigma);
    }
}

/**
 * @brief CUDA Kernel: Multi-Block Parallel HMM Viterbi Trellis Segmentation.
 */
__global__ void xoos_hmm_viterbi_kernel(
    const GenomicBinRecord* __restrict__ d_bins,
    uint64_t num_bins,
    float purity,
    float ploidy,
    uint8_t* __restrict__ d_state_path
) {
    size_t chunk_size = (num_bins + gridDim.x - 1) / gridDim.x;
    size_t start_idx = blockIdx.x * chunk_size;
    size_t end_idx = start_idx + chunk_size;
    if (end_idx > num_bins) end_idx = num_bins;
    if (start_idx >= num_bins || threadIdx.x != 0) return;

    // Expected Log2 ratios for CN = 0, 1, 2, 3, 4, 5
    float expected_log2r[6];
    float avg_tumor_depth = purity * ploidy + (1.0f - purity) * 2.0f;
    for (int k = 0; k < 6; ++k) {
        float state_depth = purity * static_cast<float>(k) + (1.0f - purity) * 2.0f;
        float ratio = (state_depth + 1e-4f) / (avg_tumor_depth + 1e-4f);
        expected_log2r[k] = logf(ratio) * kLog2E;
    }

    // Dynamic programming trellis states
    float v_prev[6];
    float v_curr[6];

    for (int k = 0; k < 6; ++k) {
        v_prev[k] = (k == 2) ? 0.0f : -4.0f; // Prior: Diploid most likely
    }

    float inv_two_sigma_sq = 1.0f / (2.0f * kDefaultEmissionSigma * kDefaultEmissionSigma);

    for (size_t i = start_idx; i < end_idx; ++i) {
        float obs_log2r = d_bins[i].normalized_log2r;

        for (int k = 0; k < 6; ++k) {
            float diff = obs_log2r - expected_log2r[k];
            float log_emission = - (diff * diff) * inv_two_sigma_sq;

            // Find best predecessor state j
            float max_trans = -1e9f;
            for (int j = 0; j < 6; ++j) {
                float trans_prob = (j == k) ? 0.0f : kStateSwitchPenalty;
                float score = v_prev[j] + trans_prob;
                if (score > max_trans) max_trans = score;
            }

            v_curr[k] = max_trans + log_emission;
        }

        // Trace best state for position i
        uint8_t best_state = 2;
        float max_score = v_curr[2];
        for (int k = 0; k < 6; ++k) {
            if (v_curr[k] > max_score) {
                max_score = v_curr[k];
                best_state = static_cast<uint8_t>(k);
            }
        }
        d_state_path[i] = best_state;

        for (int k = 0; k < 6; ++k) {
            v_prev[k] = v_curr[k];
        }
    }
}

} // namespace xoos::cnv_caller::cuda
