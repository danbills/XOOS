#pragma once

/**
 * ============================================================================
 * Short Tandem Repeat (STR) & Expansion Caller Core Types
 * File: str_caller/src/cuda/str_types.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Design Goals:
 *   - 16-byte aligned read evidence records with repeat tract lengths and alignment qualities.
 *   - Polymerase stutter noise parameters (downstream slippage d, upstream insertion u).
 *   - 2D diploid genotype matrix (Allele 1, Allele 2) for repeat lengths up to 128 units.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <cuda_runtime.h>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace xoos::str_caller::cuda {

constexpr uint32_t kMaxRepeatUnits = 128;      // Max repeat count for standard expansion scanning
constexpr uint32_t kMaxMotifLen = 16;          // Max motif length (e.g. trinucleotide CAG, hexanucleotide GGGGCC)

/**
 * @brief STR Locus Descriptor.
 */
struct __align__(16) StrLocusDescriptor {
    uint32_t locus_id;
    uint32_t chr_id;
    uint32_t start_pos;
    uint32_t end_pos;
    uint16_t motif_len;
    uint16_t ref_repeat_count;
    char motif[kMaxMotifLen];
    float down_stutter_prob;   // Slippage rate d (typically 0.05 - 0.20)
    float up_stutter_prob;     // Insertion rate u (typically 0.01 - 0.05)
    float geometric_factor;    // Decay constant for stutter distribution
};

/**
 * @brief Read Evidence supporting an STR locus.
 */
struct __align__(16) StrReadEvidence {
    uint64_t read_id;
    uint32_t locus_id;
    uint16_t observed_repeat_count;
    uint8_t read_type;         // 0=Spanning, 1=Flanking, 2=In-Repeat
    uint8_t mapq;
    float alignment_score;
    float base_qual_avg;
};

/**
 * @brief Genotype Call Result for an STR Locus.
 */
struct __align__(16) StrGenotypeCall {
    uint32_t locus_id;
    uint16_t allele1;          // Inferred shorter repeat allele
    uint16_t allele2;          // Inferred longer repeat allele
    float log_likelihood;      // Genotype log likelihood
    float genotype_quality;    // Phred-scaled GQ
    uint32_t total_support_reads;
    uint32_t spanning_reads;
    uint32_t in_repeat_reads;
    bool is_expansion;         // True if allele exceeds clinical threshold
};

/**
 * @brief STR Execution Statistics.
 */
struct StrExecutionStats {
    uint64_t total_loci_processed = 0;
    uint64_t total_reads_evaluated = 0;
    uint64_t total_genotypes_tested = 0;
    double kernel_time_ms = 0.0;
    double total_time_ms = 0.0;
    double throughput_reads_per_sec = 0.0;
};

} // namespace xoos::str_caller::cuda
