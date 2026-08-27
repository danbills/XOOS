#pragma once

/**
 * ============================================================================
 * Read Collapser Core Types & Data Layouts
 * File: read_collapser/src/cuda/collapser_types.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Design Goals:
 *   - 64-bit / 128-bit memory alignment (__align__(16)) for coalesced warp loads.
 *   - Fast in-register base representation:
 *       0=A, 1=C, 2=G, 3=T, 4=N, 5=Gap
 *   - Support for Duplex (hairpin-derived) and Simplex UMI duplicate clusters.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <cuda_runtime.h>
#include <cstdint>
#include <cstddef>

namespace xoos::read_collapser::cuda {

constexpr uint32_t kMaxClusterSize = 256;      // Maximum reads in a duplicate cluster
constexpr uint32_t kMaxCollapsedReadLen = 256;  // Maximum bases in collapsed read

enum class BaseIndex : uint8_t {
    kA = 0,
    kC = 1,
    kG = 2,
    kT = 3,
    kN = 4,
    kGap = 5,
    kInvalid = 6
};

/**
 * @brief Base character to 3-bit index conversion.
 */
__host__ __device__ __forceinline__ BaseIndex char_to_base_idx(char c) {
    switch (c) {
        case 'A': case 'a': return BaseIndex::kA;
        case 'C': case 'c': return BaseIndex::kC;
        case 'G': case 'g': return BaseIndex::kG;
        case 'T': case 't': return BaseIndex::kT;
        case 'N': case 'n': return BaseIndex::kN;
        case '-': case '.': return BaseIndex::kGap;
        default:            return BaseIndex::kN;
    }
}

/**
 * @brief 3-bit index to Base character conversion.
 */
__host__ __device__ __forceinline__ char base_idx_to_char(BaseIndex idx) {
    switch (idx) {
        case BaseIndex::kA: return 'A';
        case BaseIndex::kC: return 'C';
        case BaseIndex::kG: return 'G';
        case BaseIndex::kT: return 'T';
        case BaseIndex::kN: return 'N';
        case BaseIndex::kGap: return '-';
        default:            return 'N';
    }
}

/**
 * @brief Compact descriptor of an input aligned read for clustering.
 */
struct __align__(16) ReadClusterDescriptor {
    uint64_t umi_hash;    // 64-bit hash of UMI sequence (or SID)
    uint64_t rbeg;        // Reference alignment start coordinate
    int32_t read_idx;     // Global index in input read buffer
    int32_t read_len;     // Length of read in bases
    int32_t is_reverse;   // 1 if reverse strand, 0 if forward strand
    int32_t is_duplex;    // 1 if duplex consensus read, 0 if simplex
};

/**
 * @brief Cluster boundary interval [start_idx, end_idx).
 */
struct __align__(8) ReadClusterSpan {
    uint32_t start_idx;   // Offset in sorted ReadClusterDescriptor array
    uint32_t num_reads;   // Number of duplicate reads in this cluster
    uint64_t rbeg;        // Genomic start coordinate
    uint64_t umi_hash;    // Cluster UMI
};

/**
 * @brief Column base counts and quality accumulations for consensus solving.
 */
struct __align__(16) ColumnBaseCountsGPU {
    uint16_t fwd_counts[6];  // A, C, G, T, N, Gap
    uint16_t rev_counts[6];  // A, C, G, T, N, Gap
    uint32_t fwd_qual_sum[6]; // Phred quality sums
    uint32_t rev_qual_sum[6]; // Phred quality sums
};

/**
 * @brief Final collapsed consensus record output.
 */
struct __align__(16) CollapsedReadResult {
    uint64_t rbeg;
    uint32_t cluster_size;
    uint16_t consensus_len;
    uint16_t mean_qual;
    char consensus_seq[kMaxCollapsedReadLen];
    char consensus_qual[kMaxCollapsedReadLen];
    char yc_tag[kMaxCollapsedReadLen];
};

/**
 * @brief Batch execution statistics.
 */
struct CollapserStats {
    uint64_t total_input_reads = 0;
    uint64_t total_clusters_formed = 0;
    uint64_t total_collapsed_reads = 0;
    double clustering_time_ms = 0.0;
    double consensus_time_ms = 0.0;
    double total_time_ms = 0.0;
    double throughput_reads_per_sec = 0.0;
    double throughput_gbps = 0.0;
};

} // namespace xoos::read_collapser::cuda
