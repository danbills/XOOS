#pragma once

/**
 * ============================================================================
 * GPU Background Noise & Global Tumor Fraction Solver Kernel
 * File: tumor_fraction_estimator/src/cuda/noise_correction_solver_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Mathematical Foundations:
 *   - Global artifact rate estimation: \epsilon = \sum OtherAlts / (2 * \sum Depth).
 *   - Background noise-subtracted tumor fraction: TF = 2 * Mean(VAF_{adj}).
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include "tfe_types.cuh"
#include <cuda_runtime.h>
#include <cstdint>

namespace xoos::tfe::cuda {

/**
 * @brief CUDA Kernel: Reduce site pileups into global tumor fraction and noise summary.
 */
__global__ void xoos_tfe_summary_solver_kernel(
    const VariantSitePileupResult* __restrict__ d_results,
    uint32_t num_sites,
    TumorFractionSummary* __restrict__ d_summary
) {
    if (threadIdx.x != 0 || blockIdx.x != 0) return;

    uint64_t tot_ref = 0;
    uint64_t tot_alt = 0;
    uint64_t tot_other = 0;
    uint64_t tot_depth = 0;
    uint32_t pass_sites = 0;
    uint32_t detected = 0;

    for (uint32_t i = 0; i < num_sites; ++i) {
        const VariantSitePileupResult& r = d_results[i];
        if (r.is_passed) {
            pass_sites++;
            tot_ref += r.ref_count;
            tot_alt += r.alt_count;
            tot_other += r.other_alts_count;
            tot_depth += r.total_depth;
            if (r.alt_count > 0) detected++;
        }
    }

    double adj_alt = 0.0;
    if (tot_alt > (tot_other / 2)) {
        adj_alt = static_cast<double>(tot_alt) - 0.5 * static_cast<double>(tot_other);
    }

    double mean_vaf = (tot_depth > 0) ? (adj_alt / tot_depth) : 0.0;
    double err_rate = (tot_depth > 0) ? (static_cast<double>(tot_other) / (2.0 * tot_depth)) : 0.0;

    TumorFractionSummary sum;
    sum.tumor_fraction = 2.0 * mean_vaf;
    sum.mean_vaf = mean_vaf;
    sum.error_rate = err_rate;
    sum.total_alt = tot_alt;
    sum.total_ref = tot_ref;
    sum.total_other_alts = tot_other;
    sum.total_depth = tot_depth;
    sum.total_adjusted_alt = adj_alt;
    sum.total_sites = num_sites;
    sum.passing_sites = pass_sites;
    sum.sites_detected = detected;

    *d_summary = sum;
}

} // namespace xoos::tfe::cuda
