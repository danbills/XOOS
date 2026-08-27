#pragma once

/**
 * ============================================================================
 * Copy Number Caller CUDA Engine: Executive Header
 * File: copy_number_caller/src/cuda/cnv_caller_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Capabilities:
 *   - Fast Warp-Parallel GC Bias Normalization & Log2 Ratio Processing.
 *   - HMM Viterbi Trellis Segmentation across hundreds of thousands of bins.
 *   - 2D Grid Optimization Solver for Global Tumor Purity & Ploidy Estimation.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <vector>
#include <string>
#include <memory>
#include <cuda_runtime.h>

#include "cnv_types.cuh"

namespace xoos::cnv_caller::cuda {

class CopyNumberCallerCudaEngine {
public:
    explicit CopyNumberCallerCudaEngine(int device_id = 0);
    ~CopyNumberCallerCudaEngine();

    /**
     * @brief Perform complete GC bias correction and Log2 ratio transformation.
     */
    bool normalize_gc_and_log2r(
        std::vector<GenomicBinRecord>& bins,
        float gc_a = -0.5f,
        float gc_b = 1.0f,
        float gc_c = 50.0f,
        float baseline_depth = 100.0f
    );

    /**
     * @brief Run end-to-end CNV analysis: GC normalization, HMM segmentation, and purity/ploidy optimization.
     */
    bool run_cnv_calling(
        std::vector<GenomicBinRecord>& bins,
        std::vector<CnvSegment>& out_segments,
        PurityPloidyFit& out_fit,
        CnvExecutionStats& out_stats
    );

    /**
     * @brief Print device status and SM hardware metrics.
     */
    void print_device_info() const;

private:
    int device_id_ = 0;
    cudaStream_t stream_ = nullptr;

    GenomicBinRecord* d_bins_ = nullptr;
    uint8_t* d_state_path_ = nullptr;
    CnvSegment* d_segments_ = nullptr;
    float* d_grid_scores_ = nullptr;
    PurityPloidyFit* d_fit_result_ = nullptr;

    size_t max_bins_ = 2000000;
    size_t max_segments_ = 65536;
    size_t max_grid_points_ = kPurityGridSteps * kPloidyGridSteps;

    void allocate_workspace(size_t max_bins, size_t max_segs);
    void free_workspace();
};

} // namespace xoos::cnv_caller::cuda
