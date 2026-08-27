#pragma once

/**
 * ============================================================================
 * GPU Parallel Variant Pileup & Quality Filter Kernel
 * File: tumor_fraction_estimator/src/cuda/variant_pileup_filter_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Capabilities:
 *   - Warp-Parallel Read Quality & Base Quality Filtering.
 *   - Ref / Alt / Other-Alts Classification & Depth Aggregation.
 *   - Site-Level Pass/Fail Gate Evaluation.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include "tfe_types.cuh"
#include <cuda_runtime.h>
#include <cstdint>

namespace xoos::tfe::cuda {

/**
 * @brief CUDA Kernel: Multi-Threaded Variant Pileup & Quality Filtering.
 */
__global__ void xoos_variant_pileup_kernel(
    const VariantProbeSite* __restrict__ d_sites,
    const ProbeReadObservation* __restrict__ d_reads,
    const uint32_t* __restrict__ d_site_offsets,
    const uint32_t* __restrict__ d_site_counts,
    uint32_t num_sites,
    VariantSitePileupResult* __restrict__ d_results
) {
    uint32_t site_idx = blockIdx.x;
    if (site_idx >= num_sites) return;

    const VariantProbeSite site = d_sites[site_idx];
    uint32_t offset = d_site_offsets[site_idx];
    uint32_t count = d_site_counts[site_idx];

    // Thread-local accumulation
    uint32_t t_ref = 0;
    uint32_t t_alt = 0;
    uint32_t t_other = 0;

    for (uint32_t r = threadIdx.x; r < count; r += blockDim.x) {
        const ProbeReadObservation& read = d_reads[offset + r];

        // 1. Read-level quality filters
        if (read.is_duplicate != 0) continue;
        if (read.mapq < site.min_mapq) continue;
        if (read.base_qual < site.min_baseq) continue;

        // 2. Base classification
        char b = read.observed_base;
        if (b == site.ref_base) {
            t_ref++;
        } else if (b == site.alt_base) {
            t_alt++;
        } else if (b == 'A' || b == 'C' || b == 'G' || b == 'T') {
            t_other++;
        }
    }

    // Shared memory reduction within the block
    __shared__ uint32_t s_ref[256];
    __shared__ uint32_t s_alt[256];
    __shared__ uint32_t s_other[256];

    s_ref[threadIdx.x] = t_ref;
    s_alt[threadIdx.x] = t_alt;
    s_other[threadIdx.x] = t_other;
    __syncthreads();

    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            s_ref[threadIdx.x] += s_ref[threadIdx.x + s];
            s_alt[threadIdx.x] += s_alt[threadIdx.x + s];
            s_other[threadIdx.x] += s_other[threadIdx.x + s];
        }
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        uint32_t total_ref = s_ref[0];
        uint32_t total_alt = s_alt[0];
        uint32_t total_other = s_other[0];
        uint32_t depth = total_ref + total_alt + total_other;

        VariantSitePileupResult res;
        res.site_id = site.site_id;
        res.ref_count = total_ref;
        res.alt_count = total_alt;
        res.other_alts_count = total_other;
        res.total_depth = depth;

        if (depth >= 10) {
            res.observed_vaf = static_cast<float>(total_alt) / depth;
            float adj_alt = static_cast<float>(total_alt) - 0.5f * static_cast<float>(total_other);
            if (adj_alt < 0.0f) adj_alt = 0.0f;
            res.adjusted_vaf = adj_alt / depth;
            res.is_passed = 1;
        } else {
            res.observed_vaf = 0.0f;
            res.adjusted_vaf = 0.0f;
            res.is_passed = 0;
        }

        d_results[site_idx] = res;
    }
}

} // namespace xoos::tfe::cuda
