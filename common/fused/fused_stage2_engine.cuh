#pragma once

/**
 * ============================================================================
 * Fused Stage 2 Executive Engine: AOT & JIT Hybrid Interface
 * File: common/fused/fused_stage2_engine.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Capabilities:
 *   - Fast AOT Native Execution for Standard WGS and cfDNA Canonical Profiles.
 *   - On-Demand NVRTC JIT Dynamic Execution for Arbitrary Parameter Combos.
 *   - Zero-Copy Pre-Allocated Persistent Pinned Staging Workspace (cudaMallocHost).
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <vector>
#include <string>
#include <memory>
#include <cuda_runtime.h>
#include <cuda.h>

#include "fused_types.cuh"
#include "fused_stage2_templates.cuh"
#include "nvrtc_jit_engine.hpp"

namespace xoos::fused::cuda {

class FusedStage2CudaEngine {
public:
    explicit FusedStage2CudaEngine(int device_id = 0, size_t initial_max_reads = 3000000);
    ~FusedStage2CudaEngine();

    /**
     * @brief Get pointer to persistent pre-allocated pinned host staging memory.
     */
    FusedReadRecord* get_pinned_host_buffer(size_t required_reads);

    /**
     * @brief Run AOT pre-compiled canonical super-kernel with vector input.
     */
    bool execute_aot_canonical(
        std::vector<FusedReadRecord>& reads,
        GlobalMetricsAccumulator& out_metrics,
        FusedExecutionStats& out_stats
    );

    /**
     * @brief Run On-Demand NVRTC JIT compiled super-kernel with vector input.
     */
    bool execute_jit_dynamic(
        const jit::DynamicPolicyConfig& config,
        std::vector<FusedReadRecord>& reads,
        GlobalMetricsAccumulator& out_metrics,
        FusedExecutionStats& out_stats
    );

    /**
     * @brief Zero-overhead JIT execution directly from persistent pinned host buffer.
     */
    bool execute_jit_dynamic_pinned(
        const jit::DynamicPolicyConfig& config,
        FusedReadRecord* h_pinned_reads,
        uint64_t num_reads,
        GlobalMetricsAccumulator& out_metrics,
        FusedExecutionStats& out_stats
    );

    /**
     * @brief Print device status and SM hardware metrics.
     */
    void print_device_info() const;

private:
    int device_id_ = 0;
    cudaStream_t stream_ = nullptr;

    FusedReadRecord* d_reads_ = nullptr;
    FusedReadRecord* h_pinned_reads_ = nullptr;
    GlobalMetricsAccumulator* d_metrics_ = nullptr;

    size_t max_reads_ = 3000000;
    std::unique_ptr<jit::NvrtcSuperKernelJit> jit_compiler_;

    void allocate_workspace(size_t max_reads);
    void free_workspace();
};

} // namespace xoos::fused::cuda
