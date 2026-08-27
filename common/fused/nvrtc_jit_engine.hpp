#pragma once

/**
 * ============================================================================
 * NVIDIA NVRTC In-Memory JIT Compilation Engine for Fused Super-Kernels
 * File: common/fused/nvrtc_jit_engine.hpp
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Features:
 *   - On-Demand C++ Policy Synthesis & Compilation via NVRTC in ~100 ms.
 *   - Disk CUBIN/PTX Caching in ~/.cache/xoos/kernels/ for 0 ms subsequent launches.
 *   - CUDA Driver API Dynamic symbol loading.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <cuda.h>

namespace xoos::fused::jit {

struct DynamicPolicyConfig {
    bool enable_rescue = true;
    bool enable_collapse = true;
    bool enable_gc_metrics = true;
    bool enable_insert_metrics = true;
    uint32_t min_family_size = 3;
    uint8_t adjusted_bq = 22;
};

class NvrtcSuperKernelJit {
public:
    NvrtcSuperKernelJit();
    ~NvrtcSuperKernelJit();

    /**
     * @brief Compile or retrieve cached CUBIN for the requested policy.
     */
    bool get_or_compile_kernel(
        const DynamicPolicyConfig& config,
        CUmodule& out_module,
        CUfunction& out_function,
        double& out_compile_time_ms
    );

    /**
     * @brief Launch kernel via dynamic driver API.
     */
    bool launch_kernel(
        CUfunction func,
        unsigned int gridDimX,
        unsigned int blockDimX,
        CUstream stream,
        void** args
    );

private:
    std::string cache_dir_;
    bool driver_initialized_ = false;

    std::string generate_source(const DynamicPolicyConfig& config) const;
    std::string compute_hash(const DynamicPolicyConfig& config) const;
    void ensure_cache_dir();
};

} // namespace xoos::fused::jit
