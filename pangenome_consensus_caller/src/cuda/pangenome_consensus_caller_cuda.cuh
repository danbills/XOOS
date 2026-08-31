#pragma once

/**
 * ============================================================================
 * Pangenome Consensus Caller CUDA Engine: Executive Header
 * File: pangenome_consensus_caller/src/cuda/pangenome_consensus_caller_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Capabilities:
 *   - Fast In-VRAM YC Tag Decoding & Mismatch Discovery.
 *   - High-Throughput Graph Alignment Delta Rescuing (> 500 Million reads/sec).
 *   - In-Place Consensus Base & Base Quality Adjustment.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <vector>
#include <string>
#include <memory>
#include <cuda_runtime.h>

#include "pangenome_types.cuh"

namespace xoos::pangenome::cuda {

class PangenomeConsensusCallerCudaEngine {
public:
    explicit PangenomeConsensusCallerCudaEngine(int device_id = 0);
    ~PangenomeConsensusCallerCudaEngine();

    /**
     * @brief Run parallel pangenome consensus rescuing across duplex reads.
     */
    bool rescue_consensus_reads(
        std::vector<DuplexReadRecord>& reads,
        std::vector<PangenomeUpdateResult>& out_results,
        PangenomeExecutionStats& out_stats
    );

    /**
     * @brief Print device status and SM hardware metrics.
     */
    void print_device_info() const;

private:
    int device_id_ = 0;
    cudaStream_t stream_ = nullptr;

    DuplexReadRecord* d_reads_ = nullptr;
    PangenomeUpdateResult* d_results_ = nullptr;

    size_t max_reads_ = 2000000;

    void allocate_workspace(size_t max_reads);
    void free_workspace();
};

} // namespace xoos::pangenome::cuda
