#pragma once

/**
 * ============================================================================
 * GPU In-Register YC Tag Lossless Duplex Decoder
 * File: pangenome_consensus_caller/src/cuda/yc_decoder_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Mathematical Foundations:
 *   - Lossless Duplex ASCII YC tag decoding in registers.
 *   - Extracts R1/R2 discrepancy alleles for graph alignment comparison.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include "pangenome_types.cuh"
#include <cuda_runtime.h>
#include <cstdint>

namespace xoos::pangenome::cuda {

/**
 * @brief Decode YC character to get alternative R2 base.
 */
__device__ __forceinline__ char decode_r2_base_from_yc(char yc, char r1_base) {
    if (yc == '*' || yc == '~' || yc == '\0') {
        return r1_base; // No discrepancy
    }

    // SBX YC substitution encoding
    // Standard mismatch codes: 'c'->(A,C), 'g'->(A,G), 't'->(A,T), 'a'->(C,A), etc.
    switch (yc) {
        case 'c': return (r1_base == 'A') ? 'C' : (r1_base == 'C' ? 'A' : r1_base);
        case 'g': return (r1_base == 'A') ? 'G' : (r1_base == 'G' ? 'A' : r1_base);
        case 't': return (r1_base == 'A') ? 'T' : (r1_base == 'T' ? 'A' : r1_base);
        case 'k': return (r1_base == 'C') ? 'G' : (r1_base == 'G' ? 'C' : r1_base);
        case 'y': return (r1_base == 'C') ? 'T' : (r1_base == 'T' ? 'C' : r1_base);
        case 'w': return (r1_base == 'G') ? 'T' : (r1_base == 'T' ? 'G' : r1_base);
        // Uppercase: inverted strand or indel variants
        case 'C': return 'C';
        case 'G': return 'G';
        case 'T': return 'T';
        case 'A': return 'A';
        default: return r1_base;
    }
}

} // namespace xoos::pangenome::cuda
