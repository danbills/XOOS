#include "aligner_cuda.cuh"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>

namespace xoos::aligner::cuda {

namespace {

/**
 * @brief CUDA Kernel: Warp-Parallel Full-Genome Read Alignment on sm_120
 */
__global__ void xoos_align_batch_kernel(
    const char* __restrict__ d_query_buffer,
    const int32_t* __restrict__ d_query_offsets,
    const int32_t* __restrict__ d_query_lens,
    uint32_t num_reads,
    const GpuFmIndex* __restrict__ d_idx,
    const char* __restrict__ d_ref_seq,
    mem_alnreg_t_GPU* __restrict__ d_out_alns,
    uint32_t* __restrict__ d_total_aligned,
    uint64_t* __restrict__ d_total_seeds,
    uint64_t* __restrict__ d_total_chains
) {
    size_t tid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    size_t stride = static_cast<size_t>(gridDim.x) * blockDim.x;

    for (size_t r = tid; r < num_reads; r += stride) {
        int32_t offset = d_query_offsets[r];
        int32_t len = d_query_lens[r];
        const char* read_seq = d_query_buffer + offset;

        // 1. Generate exact match seeds using BWT backward search
        mem_seed_t_GPU seeds[32];
        int32_t num_seeds = 0;
        generate_read_seeds_gpu(*d_idx, read_seq, len, 19, seeds, num_seeds);

        // 2. 1D Dynamic Programming DAG seed chaining
        mem_chain_t_GPU chains[8];
        int32_t num_chains = 0;
        chain_seeds_gpu(seeds, num_seeds, chains, num_chains);

        // 3. Banded DP and CIGAR generation
        mem_alnreg_t_GPU aln;
        if (num_chains > 0) {
            align_chain_to_reference_gpu(read_seq, len, d_ref_seq, chains[0], aln);
            atomicAdd(d_total_aligned, 1U);
        } else {
            // Unmapped read
            aln.rbeg = 0;
            aln.rend = 0;
            aln.score = 0;
            aln.mapq = 0;
            aln.flag = 0x04; // Unmapped bitflag
        }

        d_out_alns[r] = aln;

        // Warp-level atomic updates
        if (num_seeds > 0) {
            atomicAdd(reinterpret_cast<unsigned long long*>(d_total_seeds), static_cast<unsigned long long>(num_seeds));
        }
        if (num_chains > 0) {
            atomicAdd(reinterpret_cast<unsigned long long*>(d_total_chains), static_cast<unsigned long long>(num_chains));
        }
    }
}

} // anonymous namespace

AlignerCudaEngine::AlignerCudaEngine(int device_id) : device_id_(device_id) {
    cudaSetDevice(device_id_);
    cudaStreamCreate(&stream_);
    cudaMalloc(&d_fm_index_, sizeof(GpuFmIndex));
    allocate_workspace(max_batch_reads_, 256 * 1024 * 1024);
}

AlignerCudaEngine::~AlignerCudaEngine() {
    free_workspace();
    if (d_fm_index_) { cudaFree(d_fm_index_); d_fm_index_ = nullptr; }
    if (d_ref_seq_) { cudaFree(d_ref_seq_); d_ref_seq_ = nullptr; }
    if (stream_) { cudaStreamDestroy(stream_); stream_ = nullptr; }
}

void AlignerCudaEngine::allocate_workspace(size_t max_reads, size_t buffer_size) {
    max_batch_reads_ = max_reads;
    d_query_buffer_capacity_ = buffer_size;

    cudaMalloc(&d_query_buffer_, d_query_buffer_capacity_);
    cudaMalloc(&d_query_offsets_, max_batch_reads_ * sizeof(int32_t));
    cudaMalloc(&d_query_lens_, max_batch_reads_ * sizeof(int32_t));
    cudaMalloc(&d_aln_results_, max_batch_reads_ * sizeof(mem_alnreg_t_GPU));
}

void AlignerCudaEngine::free_workspace() {
    if (d_query_buffer_) { cudaFree(d_query_buffer_); d_query_buffer_ = nullptr; }
    if (d_query_offsets_) { cudaFree(d_query_offsets_); d_query_offsets_ = nullptr; }
    if (d_query_lens_) { cudaFree(d_query_lens_); d_query_lens_ = nullptr; }
    if (d_aln_results_) { cudaFree(d_aln_results_); d_aln_results_ = nullptr; }
}

bool AlignerCudaEngine::load_reference(
    const std::string& ref_name,
    const std::string& ref_seq,
    uint64_t bwt_len,
    const std::vector<BwtOccBlock>& occ_blocks,
    const std::vector<uint64_t>& sa_table
) {
    ref_seq_len_ = ref_seq.length();
    cudaMalloc(&d_ref_seq_, ref_seq_len_);
    cudaMemcpyAsync(d_ref_seq_, ref_seq.data(), ref_seq_len_, cudaMemcpyHostToDevice, stream_);

    h_fm_index_.bwt_len = bwt_len;
    h_fm_index_.primary = 0;
    
    // Allocate device Occ blocks and SA tables
    BwtOccBlock* d_occ = nullptr;
    uint64_t* d_sa = nullptr;

    if (!occ_blocks.empty()) {
        cudaMalloc(&d_occ, occ_blocks.size() * sizeof(BwtOccBlock));
        cudaMemcpyAsync(d_occ, occ_blocks.data(), occ_blocks.size() * sizeof(BwtOccBlock), cudaMemcpyHostToDevice, stream_);
    }
    if (!sa_table.empty()) {
        cudaMalloc(&d_sa, sa_table.size() * sizeof(uint64_t));
        cudaMemcpyAsync(d_sa, sa_table.data(), sa_table.size() * sizeof(uint64_t), cudaMemcpyHostToDevice, stream_);
    }

    h_fm_index_.d_occ_table = d_occ;
    h_fm_index_.d_sa = d_sa;

    // Cumulative count initialization
    h_fm_index_.C[0] = 0;
    h_fm_index_.C[1] = 0;
    h_fm_index_.C[2] = bwt_len / 4;
    h_fm_index_.C[3] = bwt_len / 2;
    h_fm_index_.C[4] = (3 * bwt_len) / 4;

    cudaMemcpyAsync(d_fm_index_, &h_fm_index_, sizeof(GpuFmIndex), cudaMemcpyHostToDevice, stream_);
    cudaStreamSynchronize(stream_);

    std::cout << "[AlignerCudaEngine] Reference '" << ref_name << "' (" 
              << (ref_seq_len_ / 1e6) << " Mb) loaded and pinned to VRAM." << std::endl;
    return true;
}

bool AlignerCudaEngine::align_reads(
    const std::vector<std::string>& h_raw_reads,
    AlignerBatchStats& out_stats
) {
    uint32_t num_reads = static_cast<uint32_t>(h_raw_reads.size());
    if (num_reads == 0) return false;

    // Prepare contiguous read buffer
    std::string flat_buffer;
    std::vector<int32_t> offsets(num_reads);
    std::vector<int32_t> lengths(num_reads);

    size_t cur_offset = 0;
    for (uint32_t i = 0; i < num_reads; ++i) {
        offsets[i] = static_cast<int32_t>(cur_offset);
        lengths[i] = static_cast<int32_t>(h_raw_reads[i].length());
        flat_buffer += h_raw_reads[i];
        cur_offset += h_raw_reads[i].length();
    }

    if (flat_buffer.size() > d_query_buffer_capacity_) {
        free_workspace();
        allocate_workspace(num_reads * 2, flat_buffer.size() * 2);
    }

    uint32_t* d_total_aligned = nullptr;
    uint64_t* d_total_seeds = nullptr;
    uint64_t* d_total_chains = nullptr;

    cudaMalloc(&d_total_aligned, sizeof(uint32_t));
    cudaMalloc(&d_total_seeds, sizeof(uint64_t));
    cudaMalloc(&d_total_chains, sizeof(uint64_t));

    cudaMemsetAsync(d_total_aligned, 0, sizeof(uint32_t), stream_);
    cudaMemsetAsync(d_total_seeds, 0, sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_total_chains, 0, sizeof(uint64_t), stream_);

    auto t_start = std::chrono::high_resolution_clock::now();

    // DMA upload reads
    cudaMemcpyAsync(d_query_buffer_, flat_buffer.data(), flat_buffer.size(), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_query_offsets_, offsets.data(), num_reads * sizeof(int32_t), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_query_lens_, lengths.data(), num_reads * sizeof(int32_t), cudaMemcpyHostToDevice, stream_);

    // Launch alignment kernel across 170 SMs
    int threads_per_block = 256;
    int num_blocks = 170 * 4;

    xoos_align_batch_kernel<<<num_blocks, threads_per_block, 0, stream_>>>(
        d_query_buffer_, d_query_offsets_, d_query_lens_, num_reads,
        d_fm_index_, d_ref_seq_, d_aln_results_,
        d_total_aligned, d_total_seeds, d_total_chains
    );

    cudaStreamSynchronize(stream_);
    auto t_end = std::chrono::high_resolution_clock::now();

    uint32_t h_total_aligned = 0;
    uint64_t h_total_seeds = 0;
    uint64_t h_total_chains = 0;

    cudaMemcpy(&h_total_aligned, d_total_aligned, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_total_seeds, d_total_seeds, sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_total_chains, d_total_chains, sizeof(uint64_t), cudaMemcpyDeviceToHost);

    cudaFree(d_total_aligned);
    cudaFree(d_total_seeds);
    cudaFree(d_total_chains);

    double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    out_stats.total_time_ms = elapsed_ms;
    out_stats.total_reads_aligned = h_total_aligned;
    out_stats.total_seeds_generated = h_total_seeds;
    out_stats.total_chains_formed = h_total_chains;

    if (elapsed_ms > 0) {
        out_stats.throughput_reads_per_sec = num_reads / (elapsed_ms / 1000.0);
        out_stats.throughput_gbps = (flat_buffer.size() / 1e9) / (elapsed_ms / 1000.0);
    }

    return true;
}

void AlignerCudaEngine::print_device_info() const {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device_id_);
    std::cout << "\n================================================================================" << std::endl;
    std::cout << "  XOOS ALIGNER CUDA ACCELERATION ENGINE (sm_" << prop.major << prop.minor << ")" << std::endl;
    std::cout << "================================================================================" << std::endl;
    std::cout << "  Device Name:             " << prop.name << std::endl;
    std::cout << "  Compute Capability:      " << prop.major << "." << prop.minor << std::endl;
    std::cout << "  Streaming Multiprocessors: " << prop.multiProcessorCount << " SMs" << std::endl;
    std::cout << "  Total Global Memory:     " << (prop.totalGlobalMem >> 20) << " MB" << std::endl;
    std::cout << "  L2 Cache Size:           " << (prop.l2CacheSize >> 20) << " MB" << std::endl;
    std::cout << "================================================================================\n" << std::endl;
}

} // namespace xoos::aligner::cuda
