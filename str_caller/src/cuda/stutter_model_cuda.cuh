#pragma once

/**
 * ============================================================================
 * GPU Polymerase Stutter & Read Likelihood Model
 * File: str_caller/src/cuda/stutter_model_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Mathematical Foundations:
 *   - Asymmetric geometric stutter noise model for PCR/sequencing slippage.
 *   - Spanning, flanking, and in-repeat read likelihood estimators.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include "str_types.cuh"
#include <cuda_runtime.h>
#include <cmath>

namespace xoos::str_caller::cuda {

/**
 * @brief Calculate probability of observing repeat length k given true allele A under stutter noise.
 */
__device__ __forceinline__ float calculate_stutter_prob(
    uint16_t observed_k,
    uint16_t true_A,
    float down_prob,
    float up_prob,
    float rho
) {
    if (observed_k == true_A) {
        float exact_prob = 1.0f - down_prob - up_prob;
        return exact_prob > 1e-4f ? exact_prob : 1e-4f;
    } else if (observed_k < true_A) {
        int diff = true_A - observed_k;
        float decay = powf(rho, static_cast<float>(diff - 1));
        return down_prob * (1.0f - rho) * decay + 1e-6f;
    } else {
        int diff = observed_k - true_A;
        float decay = powf(rho, static_cast<float>(diff - 1));
        return up_prob * (1.0f - rho) * decay + 1e-6f;
    }
}

/**
 * @brief Calculate single-read likelihood P(Read | True Allele A).
 */
__device__ __forceinline__ float evaluate_read_likelihood(
    const StrReadEvidence& read,
    uint16_t true_A,
    const StrLocusDescriptor& locus
) {
    if (read.read_type == 0) {
        // Spanning Read: direct observation subject to stutter
        return calculate_stutter_prob(
            read.observed_repeat_count,
            true_A,
            locus.down_stutter_prob,
            locus.up_stutter_prob,
            locus.geometric_factor
        );
    } else if (read.read_type == 1) {
        // Flanking Read: partial observation (at least k repeats)
        if (true_A >= read.observed_repeat_count) {
            return 0.75f;
        } else {
            return 0.05f;
        }
    } else {
        // In-Repeat Read: long expansion support (A must exceed read capacity)
        if (true_A >= read.observed_repeat_count) {
            return 0.90f;
        } else {
            return 0.01f;
        }
    }
}

} // namespace xoos::str_caller::cuda
