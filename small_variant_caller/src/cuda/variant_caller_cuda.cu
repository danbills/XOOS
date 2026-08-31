#include "variant_caller_cuda.cuh"
#include "pairhmm_cuda.cuh"
#include "genotype_caller_cuda.cuh"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cmath>

namespace xoos::variant_caller::cuda {

SmallVariantCallerCudaEngine::SmallVariantCallerCudaEngine(int device_id) : device_id_(device_id) {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaSetDevice(device_id_ < count ? device_id_ : 0);
    }
    cudaStreamCreate(&stream_);
    init_pairhmm_constant_tables();
    allocate_workspace(max_batch_pairs_, max_reads_, max_haps_);
}

SmallVariantCallerCudaEngine::~SmallVariantCallerCudaEngine() {
    free_workspace();
    if (stream_) { cudaStreamDestroy(stream_); stream_ = nullptr; }
}

void SmallVariantCallerCudaEngine::allocate_workspace(size_t max_pairs, size_t max_reads, size_t max_haps) {
    free_workspace();
    max_batch_pairs_ = max_pairs;
    max_reads_ = max_reads;
    max_haps_ = max_haps;

    cudaMalloc(&d_reads_, max_reads_ * sizeof(ActiveRegionRead));
    cudaMalloc(&d_haplotypes_, max_haps_ * sizeof(HaplotypeDescriptor));
    cudaMalloc(&d_read_indices_, max_batch_pairs_ * sizeof(uint32_t));
    cudaMalloc(&d_hap_indices_, max_batch_pairs_ * sizeof(uint32_t));
    cudaMalloc(&d_log_likelihoods_, max_batch_pairs_ * sizeof(double));
    cudaMalloc(&d_variant_result_, sizeof(VariantCallResult));
}

void SmallVariantCallerCudaEngine::free_workspace() {
    if (d_reads_) { cudaFree(d_reads_); d_reads_ = nullptr; }
    if (d_haplotypes_) { cudaFree(d_haplotypes_); d_haplotypes_ = nullptr; }
    if (d_read_indices_) { cudaFree(d_read_indices_); d_read_indices_ = nullptr; }
    if (d_hap_indices_) { cudaFree(d_hap_indices_); d_hap_indices_ = nullptr; }
    if (d_log_likelihoods_) { cudaFree(d_log_likelihoods_); d_log_likelihoods_ = nullptr; }
    if (d_variant_result_) { cudaFree(d_variant_result_); d_variant_result_ = nullptr; }
}

bool SmallVariantCallerCudaEngine::compute_pairhmm(
    const std::vector<ActiveRegionRead>& h_reads,
    const std::vector<HaplotypeDescriptor>& h_haplotypes,
    const std::vector<uint32_t>& h_read_indices,
    const std::vector<uint32_t>& h_hap_indices,
    std::vector<double>& out_log_likelihoods,
    VariantExecutionStats& out_exec_stats
) {
    uint64_t num_pairs = h_read_indices.size();
    if (num_pairs == 0 || h_reads.empty() || h_haplotypes.empty()) return false;

    if (num_pairs > max_batch_pairs_ || h_reads.size() > max_reads_ || h_haplotypes.size() > max_haps_) {
        allocate_workspace(std::max(num_pairs, max_batch_pairs_ * 2),
                           std::max(h_reads.size(), max_reads_ * 2),
                           std::max(h_haplotypes.size(), max_haps_ * 2));
    }

    out_log_likelihoods.resize(num_pairs);

    // 1. Transfer reads, haplotypes, indices
    cudaMemcpyAsync(d_reads_, h_reads.data(), h_reads.size() * sizeof(ActiveRegionRead), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_haplotypes_, h_haplotypes.data(), h_haplotypes.size() * sizeof(HaplotypeDescriptor), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_read_indices_, h_read_indices.data(), num_pairs * sizeof(uint32_t), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_hap_indices_, h_hap_indices.data(), num_pairs * sizeof(uint32_t), cudaMemcpyHostToDevice, stream_);

    // 2. Launch PairHMM Kernel
    cudaEvent_t k_start, k_stop;
    cudaEventCreate(&k_start);
    cudaEventCreate(&k_stop);
    cudaEventRecord(k_start, stream_);

    int threads = 128;
    int blocks = (num_pairs + threads - 1) / threads;
    if (blocks > 170 * 16) blocks = 170 * 16;
    if (blocks == 0) blocks = 1;

    xoos_pairhmm_kernel<<<blocks, threads, 0, stream_>>>(
        d_reads_, d_haplotypes_, d_read_indices_, d_hap_indices_, num_pairs, d_log_likelihoods_
    );

    cudaEventRecord(k_stop, stream_);
    cudaEventSynchronize(k_stop);

    float kernel_time_ms = 0.0f;
    cudaEventElapsedTime(&kernel_time_ms, k_start, k_stop);
    cudaEventDestroy(k_start);
    cudaEventDestroy(k_stop);

    // 3. Download results
    cudaMemcpyAsync(out_log_likelihoods.data(), d_log_likelihoods_, num_pairs * sizeof(double), cudaMemcpyDeviceToHost, stream_);
    cudaStreamSynchronize(stream_);

    // 4. Calculate total dynamic programming cells & GCUPS
    uint64_t total_cells = 0;
    for (uint64_t i = 0; i < num_pairs; ++i) {
        uint32_t r_idx = h_read_indices[i];
        uint32_t h_idx = h_hap_indices[i];
        total_cells += static_cast<uint64_t>(h_reads[r_idx].length) * h_haplotypes[h_idx].length;
    }

    out_exec_stats.total_pairhmm_matrices = num_pairs;
    out_exec_stats.total_dp_cells_computed = total_cells;
    out_exec_stats.pairhmm_time_ms = kernel_time_ms;
    out_exec_stats.total_time_ms = kernel_time_ms;
    if (kernel_time_ms > 0) {
        out_exec_stats.gcups = (static_cast<double>(total_cells) / 1e9) / (kernel_time_ms / 1000.0);
        out_exec_stats.throughput_matrices_per_sec = num_pairs / (kernel_time_ms / 1000.0);
    }

    return true;
}

bool SmallVariantCallerCudaEngine::call_variant(
    const std::vector<ActiveRegionRead>& h_reads,
    const std::vector<HaplotypeDescriptor>& h_haplotypes,
    uint64_t variant_pos,
    char ref_base,
    char alt_base,
    VariantCallResult& out_result,
    VariantExecutionStats& out_exec_stats
) {
    uint32_t num_reads = static_cast<uint32_t>(h_reads.size());
    uint32_t num_haps = static_cast<uint32_t>(h_haplotypes.size());
    if (num_reads == 0 || num_haps < 2) return false;

    // Form dense (Read, Haplotype) pairs
    uint64_t total_pairs = static_cast<uint64_t>(num_reads) * num_haps;
    std::vector<uint32_t> r_indices(total_pairs);
    std::vector<uint32_t> h_indices(total_pairs);

    for (uint32_t r = 0; r < num_reads; ++r) {
        for (uint32_t h = 0; h < num_haps; ++h) {
            uint64_t p = static_cast<uint64_t>(r) * num_haps + h;
            r_indices[p] = r;
            h_indices[p] = h;
        }
    }

    std::vector<double> log_lks;
    if (!compute_pairhmm(h_reads, h_haplotypes, r_indices, h_indices, log_lks, out_exec_stats)) {
        return false;
    }

    // Launch Genotype Solver Kernel
    xoos_genotype_solver_kernel<<<1, 1, 0, stream_>>>(
        d_log_likelihoods_, d_reads_, num_reads, num_haps, variant_pos, ref_base, alt_base, d_variant_result_
    );
    cudaStreamSynchronize(stream_);

    cudaMemcpy(&out_result, d_variant_result_, sizeof(VariantCallResult), cudaMemcpyDeviceToHost);
    return true;
}

void SmallVariantCallerCudaEngine::print_device_info() const {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, device_id_ < count ? device_id_ : 0);
        std::cout << "\n================================================================================" << std::endl;
        std::cout << "  XOOS SMALL VARIANT CALLER CUDA ENGINE (sm_" << prop.major << prop.minor << ")" << std::endl;
        std::cout << "================================================================================" << std::endl;
        std::cout << "  Device Name:             " << prop.name << std::endl;
        std::cout << "  Compute Capability:      " << prop.major << "." << prop.minor << std::endl;
        std::cout << "  Streaming Multiprocessors: " << prop.multiProcessorCount << " SMs" << std::endl;
        std::cout << "  Total Global Memory:     " << (prop.totalGlobalMem >> 20) << " MB" << std::endl;
        std::cout << "  L2 Cache Size:           " << (prop.l2CacheSize >> 20) << " MB" << std::endl;
        std::cout << "================================================================================\n" << std::endl;
    }
}

} // namespace xoos::variant_caller::cuda
