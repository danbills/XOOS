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
    explicit FusedStage2CudaEngine(int device_id = 0);
    ~FusedStage2CudaEngine();

    /**
     * @brief Run AOT pre-compiled canonical super-kernel (Zero startup latency).
     */
    bool execute_aot_canonical(
        std::vector<FusedReadRecord>& reads,
        GlobalMetricsAccumulator& out_metrics,
        FusedExecutionStats& out_stats
    );

    /**
     * @brief Run On-Demand NVRTC JIT compiled super-kernel for custom configs.
     */
    bool execute_jit_dynamic(
        const jit::DynamicPolicyConfig& config,
        std::vector<FusedReadRecord>& reads,
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
    GlobalMetricsAccumulator* d_metrics_ = nullptr;

    size_t max_reads_ = 2000000;
    std::unique_ptr<jit::NvrtcSuperKernelJit> jit_compiler_;

    void allocate_workspace(size_t max_reads);
    void free_workspace();
};

} // namespace xoos::fused::cuda
