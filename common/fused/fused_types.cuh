#pragma once

/**
 * ============================================================================
 * Fused Super-Kernel Pipeline Core Types
 * File: common/fused/fused_types.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Design Goals:
 *   - 16-byte aligned unified read records bridging across pipeline stages.
 *   - Inline metrics accumulators for GC content, coverage, and insert sizes.
 *   - Fused execution statistics.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <cuda_runtime.h>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace xoos::fused::cuda {

constexpr uint32_t kMaxSeqLen = 256;
constexpr uint32_t kNumGcBins = 101; // GC% 0..100
constexpr uint32_t kMaxInsertSize = 1000;

/**
 * @brief Unified Read Record for Fused Execution.
 */
struct __align__(16) FusedReadRecord {
    uint64_t read_id;
    uint64_t barcode_hash;
    uint32_t chr_id;
    uint32_t pos;
    uint16_t length;
    uint16_t insert_size;
    uint8_t mapq;
    uint8_t is_reverse_strand;
    uint8_t has_yc_tag;
    uint8_t is_duplicate;
    float r1_graph_score;
    float r2_graph_score;
    char sequence[kMaxSeqLen];
    char yc_tag[kMaxSeqLen];
    uint8_t base_qual[kMaxSeqLen];
};

/**
 * @brief Global Inline Metrics Accumulator.
 */
struct __align__(16) GlobalMetricsAccumulator {
    uint64_t total_reads;
    uint64_t total_bases;
    uint64_t total_gc_bases;
    uint64_t total_rescued_reads;
    uint64_t total_base_corrections;
    uint64_t total_collapsed_families;
    uint64_t total_duplicates_marked;
    uint32_t gc_histogram[kNumGcBins];
    uint32_t insert_size_histogram[kMaxInsertSize];
};

/**
 * @brief Execution Statistics for Fused Super-Kernels.
 */
struct FusedExecutionStats {
    uint64_t total_reads_processed = 0;
    double jit_compile_time_ms = 0.0;
    double kernel_time_ms = 0.0;
    double total_time_ms = 0.0;
    double throughput_reads_per_sec = 0.0;
    double vram_bandwidth_gb_per_sec = 0.0;
    bool used_jit = false;
};

} // namespace xoos::fused::cuda
