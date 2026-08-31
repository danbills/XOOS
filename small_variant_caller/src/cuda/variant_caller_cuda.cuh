#pragma once

/**
 * ============================================================================
 * Small Variant Caller CUDA Engine: Executive Header
 * File: small_variant_caller/src/cuda/variant_caller_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Capabilities:
 *   - High-Throughput PairHMM Dynamic Programming (GCUPS acceleration).
 *   - Bayesian Diploid & Somatic Genotype Likelihood Solver.
 *   - Zero-Copy & Pinned VRAM Stream Execution.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <vector>
#include <string>
#include <memory>
#include <cuda_runtime.h>

#include "variant_types.cuh"

namespace xoos::variant_caller::cuda {

class SmallVariantCallerCudaEngine {
public:
    explicit SmallVariantCallerCudaEngine(int device_id = 0);
    ~SmallVariantCallerCudaEngine();

    /**
     * @brief Compute PairHMM dynamic programming likelihoods for batch of (Read, Haplotype) pairs.
     */
    bool compute_pairhmm(
        const std::vector<ActiveRegionRead>& h_reads,
        const std::vector<HaplotypeDescriptor>& h_haplotypes,
        const std::vector<uint32_t>& h_read_indices,
        const std::vector<uint32_t>& h_hap_indices,
        std::vector<double>& out_log_likelihoods,
        VariantExecutionStats& out_exec_stats
    );

    /**
     * @brief Call variant and evaluate Bayesian genotype for an active region.
     */
    bool call_variant(
        const std::vector<ActiveRegionRead>& h_reads,
        const std::vector<HaplotypeDescriptor>& h_haplotypes,
        uint64_t variant_pos,
        char ref_base,
        char alt_base,
        VariantCallResult& out_result,
        VariantExecutionStats& out_exec_stats
    );

    /**
     * @brief Print device status and SM hardware metrics.
     */
    void print_device_info() const;

private:
    int device_id_ = 0;
    cudaStream_t stream_ = nullptr;

    ActiveRegionRead* d_reads_ = nullptr;
    HaplotypeDescriptor* d_haplotypes_ = nullptr;
    uint32_t* d_read_indices_ = nullptr;
    uint32_t* d_hap_indices_ = nullptr;
    double* d_log_likelihoods_ = nullptr;
    VariantCallResult* d_variant_result_ = nullptr;

    size_t max_batch_pairs_ = 500000;
    size_t max_reads_ = 65536;
    size_t max_haps_ = 256;

    void allocate_workspace(size_t max_pairs, size_t max_reads, size_t max_haps);
    void free_workspace();
};

} // namespace xoos::variant_caller::cuda
