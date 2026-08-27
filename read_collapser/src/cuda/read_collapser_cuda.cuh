#pragma once

/**
 * ============================================================================
 * Read Collapser CUDA Engine: Top-Level Executive Interface
 * File: read_collapser/src/cuda/read_collapser_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Capabilities:
 *   - In-VRAM Parallel Duplicate Clustering (Genomic Coordinates + UMI Hash).
 *   - Intra-Warp Consensus Matrix Reduction & Majority Voting.
 *   - Phred Recalibration ($Q40 - Q60$) & IUPAC YC-Tag Generation.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <vector>
#include <string>
#include <memory>
#include <cuda_runtime.h>

#include "collapser_types.cuh"
#include "consensus_matrix_cuda.cuh"

namespace xoos::read_collapser::cuda {

class ReadCollapserCudaEngine {
public:
    explicit ReadCollapserCudaEngine(int device_id = 0);
    ~ReadCollapserCudaEngine();

    /**
     * @brief Collapse a batch of duplicate read clusters in VRAM.
     *
     * @param h_descriptors Array of read descriptors (coordinates, UMI, strand).
     * @param h_sequences Vector of read sequences.
     * @param h_qualities Vector of Phred quality strings.
     * @param out_results Vector receiving final collapsed consensus reads.
     * @param out_stats Struct populated with throughput and latency metrics.
     * @return true on success, false on error.
     */
    bool collapse_reads(
        const std::vector<ReadClusterDescriptor>& h_descriptors,
        const std::vector<std::string>& h_sequences,
        const std::vector<std::string>& h_qualities,
        std::vector<CollapsedReadResult>& out_results,
        CollapserStats& out_stats
    );

    /**
     * @brief Query device capabilities and print hardware status.
     */
    void print_device_info() const;

private:
    int device_id_ = 0;
    cudaStream_t stream_ = nullptr;

    // Device memory buffers
    ReadClusterDescriptor* d_descriptors_ = nullptr;
    char* d_seq_buffer_ = nullptr;
    char* d_qual_buffer_ = nullptr;
    int32_t* d_offsets_ = nullptr;
    int32_t* d_lens_ = nullptr;

    ReadClusterSpan* d_clusters_ = nullptr;
    CollapsedReadResult* d_collapsed_out_ = nullptr;

    size_t max_batch_reads_ = 2000000;
    size_t max_buffer_bytes_ = 256 * 1024 * 1024;

    void allocate_workspace(size_t max_reads, size_t buffer_size);
    void free_workspace();
};

} // namespace xoos::read_collapser::cuda
