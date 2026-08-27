#pragma once

/**
 * ============================================================================
 * Cell-Free DNA (cfDNA) Tumor Fraction Estimator Core Types
 * File: tumor_fraction_estimator/src/cuda/tfe_types.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Design Goals:
 *   - 16-byte aligned variant probe site descriptors.
 *   - Per-read pileup observation records with base/mapping quality.
 *   - Background error rate and noise-adjusted tumor fraction results.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <cuda_runtime.h>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace xoos::tfe::cuda {

constexpr uint8_t kVarTypeSomatic = 0;
constexpr uint8_t kVarTypeNoiseProbe = 1;
constexpr uint8_t kVarTypePopulationSNP = 2;

/**
 * @brief Variant Probe Site Descriptor.
 */
struct __align__(16) VariantProbeSite {
    uint32_t site_id;
    uint32_t chr_id;
    uint32_t pos;
    char ref_base;
    char alt_base;
    uint8_t var_type;          // Somatic, NoiseProbe, PopulationSNP
    uint8_t min_mapq;
    uint8_t min_baseq;
};

/**
 * @brief Read Observation overlapping a probe site.
 */
struct __align__(16) ProbeReadObservation {
    uint64_t read_id;
    uint32_t site_id;
    char observed_base;
    uint8_t base_qual;
    uint8_t mapq;
    uint8_t is_reverse_strand;
    uint8_t is_duplicate;
};

/**
 * @brief Pileup and Allele Counts for a single site.
 */
struct __align__(16) VariantSitePileupResult {
    uint32_t site_id;
    uint32_t ref_count;
    uint32_t alt_count;
    uint32_t other_alts_count;
    uint32_t total_depth;
    float observed_vaf;
    float adjusted_vaf;
    uint8_t is_passed;
};

/**
 * @brief Global Tumor Fraction & Quality Metrics Summary.
 */
struct __align__(16) TumorFractionSummary {
    double tumor_fraction;             // 2 * mean_vaf
    double mean_vaf;                   // Mean VAF across passing somatic sites
    double error_rate;                 // Background sequencing error rate
    uint64_t total_alt;
    uint64_t total_ref;
    uint64_t total_other_alts;
    uint64_t total_depth;
    double total_adjusted_alt;
    uint32_t total_sites;
    uint32_t passing_sites;
    uint32_t sites_detected;
};

/**
 * @brief TFE Execution Statistics.
 */
struct TfeExecutionStats {
    uint64_t total_sites_evaluated = 0;
    uint64_t total_reads_processed = 0;
    double pileup_time_ms = 0.0;
    double solver_time_ms = 0.0;
    double total_time_ms = 0.0;
    double throughput_reads_per_sec = 0.0;
};

} // namespace xoos::tfe::cuda
