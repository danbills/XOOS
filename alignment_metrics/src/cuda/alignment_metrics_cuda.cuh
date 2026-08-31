#pragma once

/**
 * ============================================================================
 * Alignment Metrics CUDA Engine: Top-Level Executive Interface
 * File: alignment_metrics/src/cuda/alignment_metrics_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Capabilities:
 *   - In-VRAM Whole-Genome / Target Depth-of-Coverage Pileup & Histogram.
 *   - Warp-Parallel Homopolymer Length Error Stratification ($1 \le L \le 16$).
 *   - Full Read Mapping Statistics & Uniformity Computation.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <vector>
#include <string>
#include <memory>
#include <cuda_runtime.h>

#include "metrics_types.cuh"

namespace xoos::alignment_metrics::cuda {

class AlignmentMetricsCudaEngine {
public:
    explicit AlignmentMetricsCudaEngine(int device_id = 0);
    ~AlignmentMetricsCudaEngine();

    /**
     * @brief Upload reference genome into GPU VRAM (pinned to L2 cache).
     */
    bool load_reference_genome(const std::string& ref_sequence);

    /**
     * @brief Compute complete alignment metrics suite in GPU VRAM.
     */
    bool compute_metrics(
        const std::vector<AlignedReadRecord>& h_reads,
        const std::vector<std::string>& h_sequences,
        CoverageSummaryMetrics& out_coverage,
        HpAccuracyMetrics& out_hp,
        ReadAlignmentStats& out_stats,
        MetricsExecutionStats& out_exec_stats
    );

    /**
     * @brief Print device status and SM hardware metrics.
     */
    void print_device_info() const;

private:
    int device_id_ = 0;
    cudaStream_t stream_ = nullptr;

    char* d_ref_genome_ = nullptr;
    uint64_t ref_length_ = 0;

    AlignedReadRecord* d_reads_ = nullptr;
    char* d_seq_buffer_ = nullptr;
    int32_t* d_offsets_ = nullptr;
    int32_t* d_lens_ = nullptr;

    int32_t* d_diff_array_ = nullptr;
    uint32_t* d_depth_array_ = nullptr;
    int32_t* d_block_sums_ = nullptr;

    uint64_t* d_coverage_hist_ = nullptr;
    uint64_t* d_hp_total_ = nullptr;
    uint64_t* d_hp_ins_ = nullptr;
    uint64_t* d_hp_del_ = nullptr;
    uint64_t* d_hp_sub_ = nullptr;

    size_t max_batch_reads_ = 2000000;
    size_t max_buffer_bytes_ = 256 * 1024 * 1024;

    void allocate_workspace(size_t max_reads, size_t buffer_size, size_t max_ref_len);
    void free_workspace();
};

} // namespace xoos::alignment_metrics::cuda
