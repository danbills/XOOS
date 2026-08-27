#pragma once

/**
 * ============================================================================
 * GPU Parallel Pangenome Graph Consensus Rescue Kernel
 * File: pangenome_consensus_caller/src/cuda/pangenome_rescuer_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Mathematical Foundations:
 *   - Parallel Evaluation of Graph Alignment Delta: \Delta = Score(R2) - Score(R1).
 *   - In-Register Consensus Allele Swapping & BQ Calibration.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include "pangenome_types.cuh"
#include "yc_decoder_cuda.cuh"
#include <cuda_runtime.h>
#include <cstdint>

namespace xoos::pangenome::cuda {

/**
 * @brief CUDA Kernel: Parallel Consensus Rescuing across millions of duplex reads.
 */
__global__ void xoos_pangenome_rescue_kernel(
    DuplexReadRecord* __restrict__ d_reads,
    uint64_t num_reads,
    PangenomeUpdateResult* __restrict__ d_results
) {
    size_t idx = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    size_t stride = static_cast<size_t>(gridDim.x) * blockDim.x;

    for (size_t r = idx; r < num_reads; r += stride) {
        DuplexReadRecord& read = d_reads[r];
        PangenomeUpdateResult res;
        res.read_id = read.read_id;
        res.num_corrections = 0;
        res.is_modified = 0;
        res.final_graph_score = read.r1_graph_score;

        // If R2 base improves pangenome graph alignment score
        if (read.has_yc_tag && (read.r2_graph_score > read.r1_graph_score)) {
            uint16_t len = read.length;
            if (len > kMaxReadLen) len = kMaxReadLen;

            for (uint16_t i = 0; i < len; ++i) {
                char yc = read.yc_tag[i];
                if (yc != '*' && yc != '~' && yc != '\0') {
                    char r2_base = decode_r2_base_from_yc(yc, read.sequence[i]);
                    if (r2_base != read.sequence[i]) {
                        read.sequence[i] = r2_base;
                        read.base_qual[i] = kDefaultAdjustedBq;
                        res.num_corrections++;
                    }
                }
            }

            if (res.num_corrections > 0) {
                res.is_modified = 1;
                res.final_graph_score = read.r2_graph_score;
            }
        }

        d_results[r] = res;
    }
}

} // namespace xoos::pangenome::cuda
