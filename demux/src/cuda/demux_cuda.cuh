#pragma once

/**
 * ============================================================================
 * Demux CUDA Engine: Top-Level Host/Device Executive Interface
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Source Correspondence:
 *   - Directly mirrors: XOOS/demux/src/core/demux-and-trim-pipeline.h
 *   - Directly mirrors: XOOS/demux/src/task/demux.h
 *   - Directly mirrors: XOOS/demux/src/task/flow-manager.h
 *   - Directly mirrors: XOOS/demux/src/cli/demux.cpp
 *
 * Capabilities:
 *   - Warp-parallel batch demultiplexing over millions of reads in VRAM.
 *   - Supports both Duplex (SBX-D) and Simplex (YS/YSU) adapter architectures.
 *   - Zero CPU intermediate disk serialization.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <vector>
#include <string>
#include <memory>
#include <cuda_runtime.h>

#include "bitap_cuda.cuh"
#include "hairpin_cuda.cuh"
#include "duplex_consensus_cuda.cuh"

namespace xoos::demux::cuda {

struct CudaDemuxStats {
    uint64_t total_reads_processed = 0;
    uint64_t valid_duplex_reads = 0;
    uint64_t failed_hairpin_reads = 0;
    uint64_t sample_assigned_reads = 0;
    double kernel_time_ms = 0.0;
    double wallclock_time_sec = 0.0;
    double throughput_reads_per_sec = 0.0;
    double throughput_gbps = 0.0;
    float mean_concordance_rate = 0.0f;
};

class DemuxCudaEngine {
public:
    explicit DemuxCudaEngine(int device_id = 0);
    ~DemuxCudaEngine();

    /**
     * @brief Configure adapter sequences and sample index (SID) barcodes.
     */
    bool configure_adapter_bundle(
        const std::string& loop_seq,
        const std::string& start_adapter,
        const std::string& end_adapter,
        const std::vector<std::string>& sid_5p_list,
        const std::vector<std::string>& sid_3p_list
    );

    /**
     * @brief Process an in-memory batch of raw FASTQ or RDB byte records.
     *
     * @param h_raw_buffer Host pointer to raw ASCII FASTQ text buffer.
     * @param buffer_size Size in bytes of the buffer.
     * @param out_stats Struct populated with execution metrics and counters.
     * @return true on success, false on error.
     */
    bool process_raw_buffer(
        const char* h_raw_buffer,
        size_t buffer_size,
        CudaDemuxStats& out_stats
    );

    /**
     * @brief Query device capabilities and print hardware status.
     */
    void print_device_info() const;

private:
    int device_id_ = 0;
    cudaStream_t stream_ = nullptr;
    GpuAdapterBundle h_bundle_;
    GpuAdapterBundle* d_bundle_ = nullptr;
    
    // Internal streaming device workspace
    char* d_raw_buffer_ = nullptr;
    size_t d_raw_buffer_capacity_ = 0;
    
    DuplexTrimResult* d_trim_results_ = nullptr;
    ConsensusReadResult* d_consensus_results_ = nullptr;
    size_t max_batch_reads_ = 2000000; // 2M reads per chunk
    
    void allocate_workspace(size_t buffer_capacity, size_t max_reads);
    void free_workspace();
};

} // namespace xoos::demux::cuda
