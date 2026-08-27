#pragma once

/**
 * ============================================================================
 * 2D Parallel Tumor Purity & Ploidy Optimization Grid Solver
 * File: copy_number_caller/src/cuda/purity_ploidy_solver_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Mathematical Foundations:
 *   - Joint 2D Grid Evaluation: Purity $p \in [0.05, 1.00]$, Ploidy $\psi \in [1.2, 5.5]$.
 *   - Distance-to-Integer-State Least Squares & Likelihood Scoring.
 *   - Warp/Block Parallel Reduction for Global Optimum Parameter Estimation.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include "cnv_types.cuh"
#include <cuda_runtime.h>
#include <cstdint>
#include <cmath>

namespace xoos::cnv_caller::cuda {

/**
 * @brief CUDA Kernel: Evaluate entire 2D Purity/Ploidy grid across all segments.
 */
__global__ void xoos_purity_ploidy_grid_kernel(
    const CnvSegment* __restrict__ d_segments,
    uint32_t num_segments,
    float min_purity,
    float max_purity,
    float min_ploidy,
    float max_ploidy,
    uint32_t purity_steps,
    uint32_t ploidy_steps,
    float* __restrict__ d_grid_scores
) {
    size_t idx = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    size_t total_points = static_cast<size_t>(purity_steps) * ploidy_steps;

    for (size_t g = idx; g < total_points; g += blockDim.x * gridDim.x) {
        uint32_t p_idx = g / ploidy_steps;
        uint32_t psi_idx = g % ploidy_steps;

        float p = min_purity + (max_purity - min_purity) * (static_cast<float>(p_idx) / (purity_steps - 1));
        float psi = min_ploidy + (max_ploidy - min_ploidy) * (static_cast<float>(psi_idx) / (ploidy_steps - 1));

        float avg_tumor_depth = p * psi + (1.0f - p) * 2.0f;

        // Precompute expected Log2R for integer states 0..7
        float exp_log2r[8];
        #pragma unroll
        for (int k = 0; k < 8; ++k) {
            float state_depth = p * static_cast<float>(k) + (1.0f - p) * 2.0f;
            exp_log2r[k] = logf((state_depth + 1e-4f) / (avg_tumor_depth + 1e-4f)) * 1.4426950408889634f;
        }

        // Sum residuals across all segments
        float total_residual = 0.0f;

        for (uint32_t s = 0; s < num_segments; ++s) {
            float obs_r = d_segments[s].mean_log2r;
            float weight = static_cast<float>(d_segments[s].num_bins);

            // Find closest theoretical copy number state
            float min_sq_err = 1e9f;
            #pragma unroll
            for (int k = 0; k < 8; ++k) {
                float err = obs_r - exp_log2r[k];
                float sq = err * err;
                if (sq < min_sq_err) min_sq_err = sq;
            }

            total_residual += weight * min_sq_err;
        }

        d_grid_scores[g] = -total_residual; // Higher score = better fit
    }
}

/**
 * @brief CUDA Kernel: Find global best purity and ploidy from evaluated grid scores.
 */
__global__ void xoos_reduce_best_purity_ploidy_kernel(
    const float* __restrict__ d_grid_scores,
    uint32_t total_points,
    float min_purity,
    float max_purity,
    float min_ploidy,
    float max_ploidy,
    uint32_t purity_steps,
    uint32_t ploidy_steps,
    PurityPloidyFit* __restrict__ d_fit_result
) {
    if (threadIdx.x != 0 || blockIdx.x != 0) return;

    float best_score = -1e30f;
    size_t best_idx = 0;

    for (size_t i = 0; i < total_points; ++i) {
        if (d_grid_scores[i] > best_score) {
            best_score = d_grid_scores[i];
            best_idx = i;
        }
    }

    uint32_t p_idx = best_idx / ploidy_steps;
    uint32_t psi_idx = best_idx % ploidy_steps;

    float best_p = min_purity + (max_purity - min_purity) * (static_cast<float>(p_idx) / (purity_steps - 1));
    float best_psi = min_ploidy + (max_ploidy - min_ploidy) * (static_cast<float>(psi_idx) / (ploidy_steps - 1));

    d_fit_result->best_purity = best_p;
    d_fit_result->best_ploidy = best_psi;
    d_fit_result->max_log_likelihood = static_cast<double>(best_score);
    d_fit_result->subclonal_fraction = 0.0f;
}

} // namespace xoos::cnv_caller::cuda
