#pragma once

/**
 * ============================================================================
 * Short Tandem Repeat (STR) Caller CUDA Engine: Executive Header
 * File: str_caller/src/cuda/str_caller_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Capabilities:
 *   - Warp-Parallel Asymmetric Polymerase Stutter Evaluation.
 *   - Massive In-Register 2D Diploid Genotype Grid Search across thousands of loci.
 *   - Pathogenic Repeat Expansion Detection & Phred GQ Calculation.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <vector>
#include <string>
#include <memory>
#include <cuda_runtime.h>

#include "str_types.cuh"

namespace xoos::str_caller::cuda {

class StrCallerCudaEngine {
public:
    explicit StrCallerCudaEngine(int device_id = 0);
    ~StrCallerCudaEngine();

    /**
     * @brief Run parallel STR genotyping & repeat expansion calling across loci.
     */
    bool genotype_loci(
        const std::vector<StrLocusDescriptor>& loci,
        const std::vector<StrReadEvidence>& reads,
        const std::vector<uint32_t>& locus_read_offsets,
        const std::vector<uint32_t>& locus_read_counts,
        std::vector<StrGenotypeCall>& out_calls,
        StrExecutionStats& out_stats,
        uint16_t max_repeat_search = 64
    );

    /**
     * @brief Print device status and SM hardware metrics.
     */
    void print_device_info() const;

private:
    int device_id_ = 0;
    cudaStream_t stream_ = nullptr;

    StrLocusDescriptor* d_loci_ = nullptr;
    StrReadEvidence* d_reads_ = nullptr;
    uint32_t* d_offsets_ = nullptr;
    uint32_t* d_counts_ = nullptr;
    StrGenotypeCall* d_calls_ = nullptr;

    size_t max_loci_ = 65536;
    size_t max_reads_ = 2000000;

    void allocate_workspace(size_t max_loci, size_t max_reads);
    void free_workspace();
};

} // namespace xoos::str_caller::cuda
