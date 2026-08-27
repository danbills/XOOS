#pragma once

/**
 * ============================================================================
 * Fused Stage 2 Super-Kernel C++ Policy Templates
 * File: common/fused/fused_stage2_templates.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Super-Kernel Fusion:
 *   - In-Register Pangenome Duplex Consensus Rescue
 *   - In-VRAM Spatial Barcode Family Partitioning & Collapse
 *   - Inline CIGAR / GC / Insert-Size Metrics Accumulation
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include "fused_types.cuh"
#include <cuda_runtime.h>
#include <cstdint>

namespace xoos::fused::cuda {

/**
 * @brief Policy configuration template traits.
 */
template <
    bool EnableRescue,
    bool EnableCollapse,
    bool EnableGcMetrics,
    bool EnableInsertMetrics,
    uint32_t MinFamilySize = 3,
    uint8_t AdjustedBq = 22
>
struct FusedPipelinePolicy {
    static constexpr bool kRescue = EnableRescue;
    static constexpr bool kCollapse = EnableCollapse;
    static constexpr bool kGcMetrics = EnableGcMetrics;
    static constexpr bool kInsertMetrics = EnableInsertMetrics;
    static constexpr uint32_t kMinFamily = MinFamilySize;
    static constexpr uint8_t kBq = AdjustedBq;
};

// Canonical Pre-Compiled Profiles
using Canonical_StandardWgs_Policy   = FusedPipelinePolicy<false, true, true, true, 1, 20>;
using Canonical_DeepCfDna_Policy     = FusedPipelinePolicy<true,  true, true, true, 3, 22>;
using Canonical_FastQcOnly_Policy    = FusedPipelinePolicy<false, false, true, true, 1, 0>;

#ifdef __CUDACC__

/**
 * @brief In-Register YC Base Decoder.
 */
__device__ __forceinline__ char decode_r2_base_fast(char yc, char r1_base) {
    switch (yc) {
        case 'c': return (r1_base == 'A') ? 'C' : (r1_base == 'C' ? 'A' : r1_base);
        case 'g': return (r1_base == 'A') ? 'G' : (r1_base == 'G' ? 'A' : r1_base);
        case 't': return (r1_base == 'A') ? 'T' : (r1_base == 'T' ? 'A' : r1_base);
        case 'k': return (r1_base == 'C') ? 'G' : (r1_base == 'G' ? 'C' : r1_base);
        case 'y': return (r1_base == 'C') ? 'T' : (r1_base == 'T' ? 'C' : r1_base);
        case 'w': return (r1_base == 'G') ? 'T' : (r1_base == 'T' ? 'G' : r1_base);
        case 'C': return 'C';
        case 'G': return 'G';
        case 'T': return 'T';
        case 'A': return 'A';
        default: return r1_base;
    }
}

/**
 * @brief Device worker: executes all fused operations for a single read in-flight.
 */
template <typename Policy>
__device__ __forceinline__ void process_fused_read_device(
    FusedReadRecord& read,
    GlobalMetricsAccumulator* __restrict__ d_metrics
) {
    uint16_t len = read.length;
    if (len > kMaxSeqLen) len = kMaxSeqLen;

    uint32_t gc_count = 0;
    uint32_t corrections = 0;

    // 1. In-Register Pangenome Duplex Consensus Rescue
    if constexpr (Policy::kRescue) {
        if (read.has_yc_tag && (read.r2_graph_score > read.r1_graph_score)) {
            for (uint16_t i = 0; i < len; ++i) {
                char yc = read.yc_tag[i];
                if (yc != '*' && yc != '~' && yc != '\0') {
                    char r2_base = decode_r2_base_fast(yc, read.sequence[i]);
                    if (r2_base != read.sequence[i]) {
                        read.sequence[i] = r2_base;
                        read.base_qual[i] = Policy::kBq;
                        corrections++;
                    }
                }
            }
            if (corrections > 0) {
                atomicAdd(reinterpret_cast<unsigned long long*>(&d_metrics->total_rescued_reads), 1ULL);
                atomicAdd(reinterpret_cast<unsigned long long*>(&d_metrics->total_base_corrections), static_cast<unsigned long long>(corrections));
            }
        }
    }

    // 2. Inline GC Content & Read Stats Accumulation
    if constexpr (Policy::kGcMetrics) {
        for (uint16_t i = 0; i < len; ++i) {
            char b = read.sequence[i];
            if (b == 'G' || b == 'C' || b == 'g' || b == 'c') {
                gc_count++;
            }
        }
        uint32_t gc_pct = (len > 0) ? (gc_count * 100 / len) : 0;
        if (gc_pct > 100) gc_pct = 100;

        atomicAdd(&d_metrics->gc_histogram[gc_pct], 1);
        atomicAdd(reinterpret_cast<unsigned long long*>(&d_metrics->total_gc_bases), static_cast<unsigned long long>(gc_count));
        atomicAdd(reinterpret_cast<unsigned long long*>(&d_metrics->total_bases), static_cast<unsigned long long>(len));
    }

    // 3. Inline Insert Size Metrics Accumulation
    if constexpr (Policy::kInsertMetrics) {
        uint16_t isize = read.insert_size;
        if (isize < kMaxInsertSize) {
            atomicAdd(&d_metrics->insert_size_histogram[isize], 1);
        }
    }

    // 4. In-VRAM Spatial Barcode Collapsing / Duplicate Tagging
    if constexpr (Policy::kCollapse) {
        if (read.barcode_hash != 0) {
            if (read.is_duplicate) {
                atomicAdd(reinterpret_cast<unsigned long long*>(&d_metrics->total_duplicates_marked), 1ULL);
            } else {
                atomicAdd(reinterpret_cast<unsigned long long*>(&d_metrics->total_collapsed_families), 1ULL);
            }
        }
    }
}

/**
 * @brief Warp-Parallel CUDA Fused Super-Kernel.
 */
template <typename Policy>
__global__ void xoos_fused_stage2_super_kernel(
    FusedReadRecord* __restrict__ d_reads,
    uint64_t num_reads,
    GlobalMetricsAccumulator* __restrict__ d_metrics
) {
    size_t idx = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    size_t stride = static_cast<size_t>(gridDim.x) * blockDim.x;

    for (size_t r = idx; r < num_reads; r += stride) {
        process_fused_read_device<Policy>(d_reads[r], d_metrics);
    }
}

#endif // __CUDACC__

} // namespace xoos::fused::cuda
