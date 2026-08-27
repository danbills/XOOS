#pragma once

/**
 * ============================================================================
 * Aligner Core Types & Memory Layout
 * File: aligner/src/cuda/aligner_types.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Memory Alignment Rules:
 *   - Structs are packed with __align__(8) or __align__(16) to ensure 64-bit / 128-bit
 *     coalesced DRAM accesses across 32-lane warps.
 *   - No partial-sector unaligned padding to prevent read-modify-write transactions.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <cuda_runtime.h>
#include <cstdint>
#include <cstddef>

namespace xoos::aligner::cuda {

using bwtint_t = uint64_t;

/**
 * @brief BWT Search Interval on GPU (Occurrence Interval)
 */
struct __align__(16) bwtintv_t_GPU {
    bwtint_t x[3]; // x[0]=k (lower bound), x[1]=l (upper bound), x[2]=s (interval size)
    bwtint_t info; // Auxiliary metadata (query position, seed length)
};

/**
 * @brief Exact Seed Match on Reference Genome
 */
struct __align__(8) mem_seed_t_GPU {
    uint64_t rbeg;    // 0-based reference start coordinate (forward strand)
    int32_t qbeg;     // 0-based query read start coordinate
    int32_t len;      // Seed length in bases
    int32_t score;    // Seeding score
};

/**
 * @brief Chained Alignment Region (Seed Chain)
 */
struct __align__(8) mem_chain_t_GPU {
    uint64_t rbeg;    // Reference coordinate of the chain anchor
    int32_t rid;      // Reference contig/chromosome ID (0..N-1)
    int32_t score;    // Aggregated chaining score
    int32_t qbeg;     // Query start coordinate
    int32_t qend;     // Query end coordinate
    int32_t n_seeds;  // Number of collinear seeds merged in this chain
};

/**
 * @brief Full Aligned Read Record (SAM / BAM Output Format)
 */
struct __align__(16) mem_alnreg_t_GPU {
    uint64_t rbeg;        // Reference start position
    uint64_t rend;        // Reference end position
    int32_t qbeg;         // Query start position
    int32_t qend;         // Query end position
    int32_t rid;          // Reference contig index
    int32_t score;        // Smith-Waterman alignment score
    uint8_t mapq;         // Mapping quality (0..60)
    uint8_t is_rev;       // 1 if aligned to reverse strand, 0 if forward
    uint16_t flag;        // SAM binary bitflag (0x00=unpaired, 0x10=reverse, etc.)
    uint32_t cigar_offset;// Offset into global CIGAR buffer
    uint16_t cigar_len;   // Number of CIGAR operations (e.g. 60M)
};

/**
 * @brief Summary statistics of GPU alignment batch execution.
 */
struct AlignerBatchStats {
    uint64_t total_reads_aligned = 0;
    uint64_t total_seeds_generated = 0;
    uint64_t total_chains_formed = 0;
    double seeding_time_ms = 0.0;
    double chaining_time_ms = 0.0;
    double dp_cigar_time_ms = 0.0;
    double total_time_ms = 0.0;
    double throughput_reads_per_sec = 0.0;
    double throughput_gbps = 0.0;
};

} // namespace xoos::aligner::cuda
