#pragma once

/**
 * ============================================================================
 * Pangenome Consensus Caller Core Types
 * File: pangenome_consensus_caller/src/cuda/pangenome_types.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Design Goals:
 *   - 16-byte aligned Duplex read structures with YC discrepancy encodings.
 *   - Graph alignment scores comparing R1 vs R2 alleles.
 *   - Consensus base correction and base quality (BQ) adjustments.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <cuda_runtime.h>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace xoos::pangenome::cuda {

constexpr uint32_t kMaxReadLen = 256;
constexpr uint8_t kDefaultAdjustedBq = 22;

/**
 * @brief Duplex Read Record with YC mismatch annotations.
 */
struct __align__(16) DuplexReadRecord {
    uint64_t read_id;
    uint16_t length;
    uint8_t mapq;
    uint8_t has_yc_tag;
    char sequence[kMaxReadLen];
    char yc_tag[kMaxReadLen];
    uint8_t base_qual[kMaxReadLen];
    float r1_graph_score;
    float r2_graph_score;
};

/**
 * @brief Result of Pangenome Consensus Adjustment for a read.
 */
struct __align__(16) PangenomeUpdateResult {
    uint64_t read_id;
    uint16_t num_corrections;
    float final_graph_score;
    uint8_t is_modified;
};

/**
 * @brief Pangenome Execution Statistics.
 */
struct PangenomeExecutionStats {
    uint64_t total_reads_processed = 0;
    uint64_t total_corrections_made = 0;
    uint64_t total_discrepancies_evaluated = 0;
    double kernel_time_ms = 0.0;
    double total_time_ms = 0.0;
    double throughput_reads_per_sec = 0.0;
};

} // namespace xoos::pangenome::cuda
