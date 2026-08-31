#pragma once

/**
 * ============================================================================
 * Aligner CUDA Engine: Top-Level Host/Device Executive Interface
 * File: aligner/src/cuda/aligner_cuda.cuh
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Capabilities:
 *   - Unified In-VRAM Seeding -> Chaining -> Banded DP Alignment.
 *   - L2 Cache Pinning of reference index on RTX 5090 (96 MB persisting window).
 *   - Direct interface for fused Demux -> Alignment pipeline.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <vector>
#include <string>
#include <memory>
#include <cuda_runtime.h>

#include "aligner_types.cuh"
#include "fmindex_cuda.cuh"
#include "chain_cuda.cuh"
#include "banded_align_cuda.cuh"

namespace xoos::aligner::cuda {

class AlignerCudaEngine {
public:
    explicit AlignerCudaEngine(int device_id = 0);
    ~AlignerCudaEngine();

    /**
     * @brief Load reference FASTA sequence and BWT occurrence tables into VRAM.
     */
    bool load_reference(
        const std::string& ref_name,
        const std::string& ref_seq,
        uint64_t bwt_len,
        const std::vector<BwtOccBlock>& occ_blocks,
        const std::vector<uint64_t>& sa_table
    );

    /**
     * @brief Align a batch of query reads directly in VRAM.
     *
     * @param h_raw_reads Vector of read sequences (or contiguous text buffer).
     * @param out_stats Struct populated with throughput and latency metrics.
     * @return true on success, false on error.
     */
    bool align_reads(
        const std::vector<std::string>& h_raw_reads,
        AlignerBatchStats& out_stats
    );

    /**
     * @brief Query device capabilities and print hardware status.
     */
    void print_device_info() const;

private:
    int device_id_ = 0;
    cudaStream_t stream_ = nullptr;
    
    // Resident reference index
    GpuFmIndex h_fm_index_;
    GpuFmIndex* d_fm_index_ = nullptr;
    char* d_ref_seq_ = nullptr;
    size_t ref_seq_len_ = 0;
    
    // Batch workspace
    char* d_query_buffer_ = nullptr;
    size_t d_query_buffer_capacity_ = 0;
    int32_t* d_query_offsets_ = nullptr;
    int32_t* d_query_lens_ = nullptr;
    
    mem_alnreg_t_GPU* d_aln_results_ = nullptr;
    size_t max_batch_reads_ = 2000000;
    
    void allocate_workspace(size_t max_reads, size_t buffer_size);
    void free_workspace();
};

} // namespace xoos::aligner::cuda
