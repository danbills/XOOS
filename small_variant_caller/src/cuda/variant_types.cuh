#pragma once

/**
 * ============================================================================
 * Small Variant Caller & PairHMM Core Types
 * File: small_variant_caller/src/cuda/variant_types.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Design Goals:
 *   - In-register PairHMM dynamic programming state representations.
 *   - 16-byte aligned read and haplotype structures for coalesced memory access.
 *   - Bayesian genotype likelihood & somatic variant feature descriptors.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <cuda_runtime.h>
#include <cstdint>
#include <cstddef>

namespace xoos::variant_caller::cuda {

constexpr uint32_t kMaxReadLen = 256;
constexpr uint32_t kMaxHapLen = 512;
constexpr uint32_t kMaxHaplotypesPerRegion = 32;

/**
 * @brief Candidate Haplotype Descriptor.
 */
struct __align__(16) HaplotypeDescriptor {
    uint32_t hap_id;
    uint32_t length;
    uint64_t start_pos;
    char sequence[kMaxHapLen];
};

/**
 * @brief Active Region Read Record for PairHMM Evaluation.
 */
struct __align__(16) ActiveRegionRead {
    uint32_t read_id;
    uint32_t length;
    uint8_t mapq;
    uint8_t is_reverse;
    char sequence[kMaxReadLen];
    uint8_t base_qual[kMaxReadLen];
    uint8_t ins_qual[kMaxReadLen];
    uint8_t del_qual[kMaxReadLen];
    uint8_t gcp_qual[kMaxReadLen]; // Gap continuation penalty
};

/**
 * @brief PairHMM Transition Parameters (Constant across read or position).
 */
struct __align__(16) PairHmmTransitionProbabilities {
    float c_mm; // Match to Match
    float c_mi; // Match to Insertion
    float c_md; // Match to Deletion
    float c_im; // Insertion to Match
    float c_ii; // Insertion to Insertion
    float c_dm; // Deletion to Match
    float c_dd; // Deletion to Deletion
};

/**
 * @brief Output Log-Likelihood of Read r given Haplotype h.
 */
struct __align__(8) ReadHaplotypeLikelihood {
    uint32_t read_id;
    uint32_t hap_id;
    double log_likelihood;
};

/**
 * @brief Genotype Call & Variant Features.
 */
struct __align__(16) VariantCallResult {
    uint64_t pos;
    char ref_allele;
    char alt_allele;
    uint32_t depth;
    uint32_t alt_depth;
    float vaf;              // Variant allele frequency (alt_depth / depth)
    float qual;             // Phred-scaled variant quality
    float strand_bias;      // Forward/reverse strand bias ratio
    uint8_t genotype;       // 0: 0/0 (HOM_REF), 1: 0/1 (HET), 2: 1/1 (HOM_ALT)
};

/**
 * @brief Execution Performance Statistics.
 */
struct VariantExecutionStats {
    uint64_t total_pairhmm_matrices = 0;
    uint64_t total_dp_cells_computed = 0;
    double pairhmm_time_ms = 0.0;
    double genotype_time_ms = 0.0;
    double total_time_ms = 0.0;
    double gcups = 0.0; // Giga Cell Updates Per Second
    double throughput_matrices_per_sec = 0.0;
};

} // namespace xoos::variant_caller::cuda
