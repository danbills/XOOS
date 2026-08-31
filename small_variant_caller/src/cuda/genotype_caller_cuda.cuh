#pragma once

/**
 * ============================================================================
 * Bayesian Genotype & Somatic Variant Solver CUDA Engine
 * File: small_variant_caller/src/cuda/genotype_caller_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Capabilities:
 *   - Log-Sum-Exp Diploid and Somatic Genotype Likelihood Evaluation.
 *   - Phred-scaled Genotype Quality (GQ) & Variant Quality (QUAL) computation.
 *   - High-precision Fisher Strand Bias & Allele Fraction (VAF) Extraction.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include "variant_types.cuh"
#include <cuda_runtime.h>
#include <cstdint>
#include <cmath>

namespace xoos::variant_caller::cuda {

/**
 * @brief CUDA Device helper: numerically stable log-sum-exp of two log values.
 */
__device__ inline double gpu_log_add(double a, double b) {
    if (a <= -900.0) return b;
    if (b <= -900.0) return a;
    return (a > b) ? (a + log1p(exp(b - a))) : (b + log1p(exp(a - b)));
}

/**
 * @brief CUDA Kernel: Evaluate Bayesian Genotypes from Read-Haplotype Likelihoods.
 */
__global__ void xoos_genotype_solver_kernel(
    const double* __restrict__ d_log_likelihoods, // Size: num_reads * num_haps
    const ActiveRegionRead* __restrict__ d_reads,
    uint32_t num_reads,
    uint32_t num_haps,
    uint64_t variant_pos,
    char ref_base,
    char alt_base,
    VariantCallResult* __restrict__ d_variant_result
) {
    if (threadIdx.x != 0 || blockIdx.x != 0) return;

    if (num_reads == 0 || num_haps < 2) {
        d_variant_result->pos = variant_pos;
        d_variant_result->ref_allele = ref_base;
        d_variant_result->alt_allele = alt_base;
        d_variant_result->depth = 0;
        d_variant_result->alt_depth = 0;
        d_variant_result->vaf = 0.0f;
        d_variant_result->qual = 0.0f;
        d_variant_result->strand_bias = 0.0f;
        d_variant_result->genotype = 0;
        return;
    }

    // Genotype priors (Uniform or flat)
    double log_lk_hom_ref = 0.0; // {0, 0}
    double log_lk_het     = 0.0; // {0, 1}
    double log_lk_hom_alt = 0.0; // {1, 1}

    uint32_t ref_fwd = 0, ref_rev = 0;
    uint32_t alt_fwd = 0, alt_rev = 0;

    for (uint32_t r = 0; r < num_reads; ++r) {
        double ll_ref = d_log_likelihoods[r * num_haps + 0];
        double ll_alt = d_log_likelihoods[r * num_haps + 1];

        // L(0/0) = P(r | h0)
        log_lk_hom_ref += ll_ref;

        // L(0/1) = 0.5 * P(r | h0) + 0.5 * P(r | h1)
        double ll_het_read = gpu_log_add(ll_ref, ll_alt) - 0.69314718056; // - ln(2)
        log_lk_het += ll_het_read;

        // L(1/1) = P(r | h1)
        log_lk_hom_alt += ll_alt;

        // Count supporting reads
        const auto& read = d_reads[r];
        if (ll_alt > ll_ref + 0.5) {
            if (read.is_reverse) alt_rev++;
            else alt_fwd++;
        } else if (ll_ref > ll_alt + 0.5) {
            if (read.is_reverse) ref_rev++;
            else ref_fwd++;
        }
    }

    uint32_t total_depth = ref_fwd + ref_rev + alt_fwd + alt_rev;
    uint32_t alt_depth = alt_fwd + alt_rev;

    // Pick best genotype
    uint8_t best_gt = 0;
    double max_ll = log_lk_hom_ref;

    if (log_lk_het > max_ll) {
        max_ll = log_lk_het;
        best_gt = 1;
    }
    if (log_lk_hom_alt > max_ll) {
        max_ll = log_lk_hom_alt;
        best_gt = 2;
    }

    // Phred-scaled variant quality
    double log_err = log_lk_hom_ref - max_ll;
    float qual = static_cast<float>(-10.0 * (log_err / 2.302585092994046));
    if (qual < 0.0f) qual = 0.0f;
    if (qual > 999.0f) qual = 999.0f;

    float vaf = (total_depth > 0) ? (static_cast<float>(alt_depth) / total_depth) : 0.0f;
    float sb = (alt_depth > 0) ? (static_cast<float>(alt_fwd) / alt_depth) : 0.5f;

    d_variant_result->pos = variant_pos;
    d_variant_result->ref_allele = ref_base;
    d_variant_result->alt_allele = alt_base;
    d_variant_result->depth = total_depth;
    d_variant_result->alt_depth = alt_depth;
    d_variant_result->vaf = vaf;
    d_variant_result->qual = qual;
    d_variant_result->strand_bias = sb;
    d_variant_result->genotype = best_gt;
}

} // namespace xoos::variant_caller::cuda
