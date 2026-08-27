#pragma once

/**
 * ============================================================================
 * High-Throughput GPU PairHMM Dynamic Programming Engine
 * File: small_variant_caller/src/cuda/pairhmm_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Mathematical Foundations:
 *   - 3-State Hidden Markov Model (Match, Insertion, Deletion).
 *   - In-Register 2-Row Rotating Buffer for O(H) local memory footprint.
 *   - Phred-to-Probability Lookup Tables stored in GPU constant memory.
 *   - Dynamic floating point scaling preventing numerical underflow.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include "variant_types.cuh"
#include <cuda_runtime.h>
#include <cstdint>
#include <cmath>

namespace xoos::variant_caller::cuda {

__constant__ float c_phred_to_prob[128];
__constant__ float c_phred_to_error[128];

/**
 * @brief Initialize Phred conversion tables in host code.
 */
inline void init_pairhmm_constant_tables() {
    float h_prob[128];
    float h_err[128];
    for (int q = 0; q < 128; ++q) {
        float err = std::pow(10.0f, -static_cast<float>(q) / 10.0f);
        h_err[q] = err;
        h_prob[q] = 1.0f - err;
    }
    cudaMemcpyToSymbol(c_phred_to_prob, h_prob, 128 * sizeof(float));
    cudaMemcpyToSymbol(c_phred_to_error, h_err, 128 * sizeof(float));
}

/**
 * @brief CUDA Kernel: High-Throughput Thread-Parallel PairHMM Engine.
 * Each GPU thread evaluates one full (Read, Haplotype) dynamic programming matrix.
 */
__global__ void xoos_pairhmm_kernel(
    const ActiveRegionRead* __restrict__ d_reads,
    const HaplotypeDescriptor* __restrict__ d_haplotypes,
    const uint32_t* __restrict__ d_read_indices,
    const uint32_t* __restrict__ d_hap_indices,
    uint64_t num_pairs,
    double* __restrict__ d_log_likelihoods
) {
    size_t idx = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    size_t stride = static_cast<size_t>(gridDim.x) * blockDim.x;

    // Local 2-row rotating buffers for DP states
    float M_prev[kMaxHapLen + 1];
    float I_prev[kMaxHapLen + 1];
    float D_prev[kMaxHapLen + 1];

    float M_curr[kMaxHapLen + 1];
    float I_curr[kMaxHapLen + 1];
    float D_curr[kMaxHapLen + 1];

    for (size_t p = idx; p < num_pairs; p += stride) {
        uint32_t r_idx = d_read_indices[p];
        uint32_t h_idx = d_hap_indices[p];

        const ActiveRegionRead& read = d_reads[r_idx];
        const HaplotypeDescriptor& hap = d_haplotypes[h_idx];

        uint32_t r_len = read.length;
        uint32_t h_len = hap.length;

        if (r_len == 0 || h_len == 0 || r_len > kMaxReadLen || h_len > kMaxHapLen) {
            d_log_likelihoods[p] = -1000.0;
            continue;
        }

        // Initialize boundary conditions: M[0, j] = 0, I[0, j] = 0, D[0, j] = 1/H
        float init_del = 1.0f / static_cast<float>(h_len);
        for (uint32_t j = 0; j <= h_len; ++j) {
            M_prev[j] = 0.0f;
            I_prev[j] = 0.0f;
            D_prev[j] = (j > 0) ? init_del : 0.0f;
        }

        double scale_log = 0.0;

        // DP Recurrence across read length
        for (uint32_t i = 1; i <= r_len; ++i) {
            char r_base = read.sequence[i - 1];
            uint8_t bq = read.base_qual[i - 1] & 0x7F;
            uint8_t iq = read.ins_qual[i - 1] & 0x7F;
            uint8_t dq = read.del_qual[i - 1] & 0x7F;
            uint8_t gq = read.gcp_qual[i - 1] & 0x7F;

            float p_err = c_phred_to_error[bq];
            float p_match = c_phred_to_prob[bq];
            float p_sub = p_err / 3.0f;

            float c_mi = c_phred_to_error[iq];
            float c_md = c_phred_to_error[dq];
            float c_mm = 1.0f - (c_mi + c_md);

            float c_ii = c_phred_to_error[gq];
            float c_im = 1.0f - c_ii;

            float c_dd = c_phred_to_error[gq];
            float c_dm = 1.0f - c_dd;

            M_curr[0] = 0.0f;
            I_curr[0] = 0.0f;
            D_curr[0] = 0.0f;

            float row_sum = 0.0f;

            for (uint32_t j = 1; j <= h_len; ++j) {
                char h_base = hap.sequence[j - 1];
                float prior = (r_base == h_base) ? p_match : p_sub;

                // Match State: M[i, j] = prior * (M[i-1,j-1]*c_mm + I[i-1,j-1]*c_im + D[i-1,j-1]*c_dm)
                float m_val = prior * (M_prev[j - 1] * c_mm + I_prev[j - 1] * c_im + D_prev[j - 1] * c_dm);

                // Insert State: I[i, j] = M[i-1, j]*c_mi + I[i-1, j]*c_ii
                float i_val = M_prev[j] * c_mi + I_prev[j] * c_ii;

                // Delete State: D[i, j] = M[i, j-1]*c_md + D[i, j-1]*c_dd
                float d_val = M_curr[j - 1] * c_md + D_curr[j - 1] * c_dd;

                M_curr[j] = m_val;
                I_curr[j] = i_val;
                D_curr[j] = d_val;

                row_sum += (m_val + i_val + d_val);
            }

            // Normalization to prevent floating point underflow
            if (row_sum > 0.0f && row_sum < 1e-15f) {
                float inv_scale = 1.0f / row_sum;
                for (uint32_t j = 0; j <= h_len; ++j) {
                    M_curr[j] *= inv_scale;
                    I_curr[j] *= inv_scale;
                    D_curr[j] *= inv_scale;
                }
                scale_log += log(static_cast<double>(row_sum));
            }

            // Rotate rows
            for (uint32_t j = 0; j <= h_len; ++j) {
                M_prev[j] = M_curr[j];
                I_prev[j] = I_curr[j];
                D_prev[j] = D_curr[j];
            }
        }

        // Final Likelihood Sum across last row
        double total_prob = 0.0;
        for (uint32_t j = 1; j <= h_len; ++j) {
            total_prob += (M_curr[j] + I_curr[j]);
        }

        if (total_prob > 0.0) {
            d_log_likelihoods[p] = log(total_prob) + scale_log;
        } else {
            d_log_likelihoods[p] = -1000.0;
        }
    }
}

} // namespace xoos::variant_caller::cuda
