#pragma once

/**
 * ============================================================================
 * Tumor Fraction Estimator CUDA Engine: Executive Header
 * File: tumor_fraction_estimator/src/cuda/tumor_fraction_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Capabilities:
 *   - Fast Warp-Parallel Allele Pileup & Filter Engine across thousands of probe sites.
 *   - Background Error Rate Estimation & Noise-Subtracted Somatic VAF Optimization.
 *   - Cell-Free DNA (cfDNA) Tumor Fraction & Contamination Quantification.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <vector>
#include <string>
#include <memory>
#include <cuda_runtime.h>

#include "tfe_types.cuh"

namespace xoos::tfe::cuda {

class TumorFractionCudaEngine {
public:
    explicit TumorFractionCudaEngine(int device_id = 0);
    ~TumorFractionCudaEngine();

    /**
     * @brief Run parallel variant pileup and tumor fraction estimation.
     */
    bool estimate_tumor_fraction(
        const std::vector<VariantProbeSite>& sites,
        const std::vector<ProbeReadObservation>& reads,
        const std::vector<uint32_t>& site_offsets,
        const std::vector<uint32_t>& site_counts,
        std::vector<VariantSitePileupResult>& out_results,
        TumorFractionSummary& out_summary,
        TfeExecutionStats& out_stats
    );

    /**
     * @brief Print device status and SM hardware metrics.
     */
    void print_device_info() const;

private:
    int device_id_ = 0;
    cudaStream_t stream_ = nullptr;

    VariantProbeSite* d_sites_ = nullptr;
    ProbeReadObservation* d_reads_ = nullptr;
    uint32_t* d_offsets_ = nullptr;
    uint32_t* d_counts_ = nullptr;
    VariantSitePileupResult* d_results_ = nullptr;
    TumorFractionSummary* d_summary_ = nullptr;

    size_t max_sites_ = 65536;
    size_t max_reads_ = 2000000;

    void allocate_workspace(size_t max_sites, size_t max_reads);
    void free_workspace();
};

} // namespace xoos::tfe::cuda
