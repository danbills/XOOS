#pragma once

/**
 * ============================================================================
 * Alignment Metrics Core Types & Data Layouts
 * File: alignment_metrics/src/cuda/metrics_types.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Design Goals:
 *   - 16-byte aligned read descriptors for coalesced warp loads.
 *   - Support for whole-genome depth of coverage, target capture coverage,
 *     homopolymer (HP) run-length accuracy, and GC bias distributions.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <cuda_runtime.h>
#include <cstdint>
#include <cstddef>

namespace xoos::alignment_metrics::cuda {

constexpr uint32_t kMaxCoverageDepth = 1024;   // Max tracked depth bin
constexpr uint32_t kMaxHomopolymerLen = 16;    // Max tracked HP length

/**
 * @brief Aligned read descriptor for parallel metric calculations.
 */
struct __align__(16) AlignedReadRecord {
    uint64_t rbeg;        // Reference start coordinate
    uint64_t rend;        // Reference end coordinate
    int32_t read_len;     // Query sequence length
    uint16_t flag;        // SAM flag (mapped, unmapped, duplicate, etc.)
    uint8_t mapq;         // Mapping quality Phred score
    uint8_t is_reverse;   // 1 if reverse strand, 0 if forward strand
};

/**
 * @brief Summary statistics of depth of coverage.
 */
struct __align__(16) CoverageSummaryMetrics {
    uint64_t total_reference_bases = 0;
    uint64_t covered_bases = 0;
    uint64_t total_aligned_bases = 0;
    double mean_coverage = 0.0;
    double median_coverage = 0.0;
    double std_coverage = 0.0;
    double pct_bases_ge_1x = 0.0;
    double pct_bases_ge_10x = 0.0;
    double pct_bases_ge_30x = 0.0;
    double pct_bases_ge_100x = 0.0;
    double pct_bases_ge_0_2x_mean = 0.0;
    double pct_bases_ge_0_5x_mean = 0.0;
    uint32_t max_coverage = 0;
    uint64_t coverage_histogram[kMaxCoverageDepth]{};
};

/**
 * @brief Homopolymer (HP) length-stratified error distribution.
 */
struct __align__(16) HpAccuracyMetrics {
    uint64_t hp_total_bases[kMaxHomopolymerLen]{};
    uint64_t hp_insertion_errors[kMaxHomopolymerLen]{};
    uint64_t hp_deletion_errors[kMaxHomopolymerLen]{};
    uint64_t hp_substitution_errors[kMaxHomopolymerLen]{};
    double hp_error_rate[kMaxHomopolymerLen]{};
};

/**
 * @brief Overall Read & Alignment Quality Summary.
 */
struct __align__(16) ReadAlignmentStats {
    uint64_t total_reads = 0;
    uint64_t mapped_reads = 0;
    uint64_t unmapped_reads = 0;
    uint64_t duplicate_reads = 0;
    uint64_t forward_strand_reads = 0;
    uint64_t reverse_strand_reads = 0;
    double mean_mapq = 0.0;
    double mean_read_length = 0.0;
};

/**
 * @brief Batch execution performance metrics.
 */
struct MetricsExecutionStats {
    uint64_t total_reads_processed = 0;
    uint64_t total_reference_bases_scanned = 0;
    double pileup_time_ms = 0.0;
    double hp_eval_time_ms = 0.0;
    double total_time_ms = 0.0;
    double throughput_reads_per_sec = 0.0;
    double throughput_gbps = 0.0;
};

} // namespace xoos::alignment_metrics::cuda
