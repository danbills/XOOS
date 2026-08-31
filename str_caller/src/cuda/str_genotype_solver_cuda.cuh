#pragma once

/**
 * ============================================================================
 * GPU Diploid Genotype Likelihood & MAP Solver Kernel
 * File: str_caller/src/cuda/str_genotype_solver_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Mathematical Foundations:
 *   - Joint Diploid Likelihood Trellis: L(A1, A2) = \sum \ln(0.5 * P(R|A1) + 0.5 * P(R|A2)).
 *   - Fast in-register 2D allele scanning for A1 <= A2 <= MaxRepeats.
 *   - Phred-scaled Genotype Quality (GQ) & Pathogenic Expansion Classification.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include "str_types.cuh"
#include "stutter_model_cuda.cuh"
#include <cuda_runtime.h>
#include <cstdint>
#include <cmath>

namespace xoos::str_caller::cuda {

constexpr float kLog10E = 0.4342944819032518f;

/**
 * @brief CUDA Kernel: Solve diploid genotypes across multiple STR loci in parallel.
 */
__global__ void xoos_str_genotype_solver_kernel(
    const StrLocusDescriptor* __restrict__ d_loci,
    const StrReadEvidence* __restrict__ d_reads,
    const uint32_t* __restrict__ d_locus_read_offsets,
    const uint32_t* __restrict__ d_locus_read_counts,
    uint32_t num_loci,
    uint16_t max_repeat_search,
    StrGenotypeCall* __restrict__ d_calls
) {
    uint32_t locus_idx = blockIdx.x;
    if (locus_idx >= num_loci) return;

    const StrLocusDescriptor locus = d_loci[locus_idx];
    uint32_t read_offset = d_locus_read_offsets[locus_idx];
    uint32_t read_count = d_locus_read_counts[locus_idx];

    // Compute total number of diploid pairs (A1 <= A2)
    uint32_t num_alleles = max_repeat_search;
    uint32_t total_pairs = (num_alleles * (num_alleles + 1)) / 2;

    float best_log_lik = -1e30f;
    float runner_up_log_lik = -1e30f;
    uint16_t best_a1 = locus.ref_repeat_count;
    uint16_t best_a2 = locus.ref_repeat_count;

    // Distribute candidate diploid pairs among threads in the block
    for (size_t pair_idx = threadIdx.x; pair_idx < total_pairs; pair_idx += blockDim.x) {
        // Map 1D index to (A1, A2) where 1 <= A1 <= A2 <= max_repeat_search
        uint32_t a1 = 1;
        uint32_t rem = pair_idx;
        while (rem >= (num_alleles - a1 + 1) && a1 < num_alleles) {
            rem -= (num_alleles - a1 + 1);
            a1++;
        }
        uint32_t a2 = a1 + rem;

        float total_log_lik = 0.0f;

        for (uint32_t r = 0; r < read_count; ++r) {
            const StrReadEvidence& read = d_reads[read_offset + r];
            float p1 = evaluate_read_likelihood(read, static_cast<uint16_t>(a1), locus);
            float p2 = evaluate_read_likelihood(read, static_cast<uint16_t>(a2), locus);
            float p_joint = 0.5f * p1 + 0.5f * p2;
            if (p_joint < 1e-12f) p_joint = 1e-12f;
            total_log_lik += logf(p_joint);
        }

        if (total_log_lik > best_log_lik) {
            runner_up_log_lik = best_log_lik;
            best_log_lik = total_log_lik;
            best_a1 = static_cast<uint16_t>(a1);
            best_a2 = static_cast<uint16_t>(a2);
        } else if (total_log_lik > runner_up_log_lik) {
            runner_up_log_lik = total_log_lik;
        }
    }

    // Warp / Block Reduction to find block-wide optimum
    // In sm_120, each block handles 1 locus
    __shared__ float s_best_lik[256];
    __shared__ float s_runner_lik[256];
    __shared__ uint16_t s_best_a1[256];
    __shared__ uint16_t s_best_a2[256];

    s_best_lik[threadIdx.x] = best_log_lik;
    s_runner_lik[threadIdx.x] = runner_up_log_lik;
    s_best_a1[threadIdx.x] = best_a1;
    s_best_a2[threadIdx.x] = best_a2;
    __syncthreads();

    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            if (s_best_lik[threadIdx.x + s] > s_best_lik[threadIdx.x]) {
                s_runner_lik[threadIdx.x] = fmaxf(s_best_lik[threadIdx.x], s_runner_lik[threadIdx.x + s]);
                s_best_lik[threadIdx.x] = s_best_lik[threadIdx.x + s];
                s_best_a1[threadIdx.x] = s_best_a1[threadIdx.x + s];
                s_best_a2[threadIdx.x] = s_best_a2[threadIdx.x + s];
            } else if (s_best_lik[threadIdx.x + s] > s_runner_lik[threadIdx.x]) {
                s_runner_lik[threadIdx.x] = s_best_lik[threadIdx.x + s];
            }
        }
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        float max_lik = s_best_lik[0];
        float second_lik = s_runner_lik[0];
        float delta = max_lik - second_lik;
        if (delta < 0.0f) delta = 0.0f;
        float gq = delta * 10.0f * kLog10E;
        if (gq > 99.0f) gq = 99.0f;

        StrGenotypeCall call;
        call.locus_id = locus.locus_id;
        call.allele1 = s_best_a1[0];
        call.allele2 = s_best_a2[0];
        call.log_likelihood = max_lik;
        call.genotype_quality = gq;
        call.total_support_reads = read_count;
        call.spanning_reads = 0;
        call.in_repeat_reads = 0;

        for (uint32_t r = 0; r < read_count; ++r) {
            if (d_reads[read_offset + r].read_type == 0) call.spanning_reads++;
            if (d_reads[read_offset + r].read_type == 2) call.in_repeat_reads++;
        }

        // Pathogenic expansion rule: allele length exceeds reference by > 50%
        call.is_expansion = (call.allele2 > locus.ref_repeat_count * 1.5f + 5);

        d_calls[locus_idx] = call;
    }
}

} // namespace xoos::str_caller::cuda
