#include "str_caller_cuda.cuh"
#include "stutter_model_cuda.cuh"
#include "str_genotype_solver_cuda.cuh"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cmath>

namespace xoos::str_caller::cuda {

StrCallerCudaEngine::StrCallerCudaEngine(int device_id) : device_id_(device_id) {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaSetDevice(device_id_ < count ? device_id_ : 0);
    }
    cudaStreamCreate(&stream_);
    allocate_workspace(max_loci_, max_reads_);
}

StrCallerCudaEngine::~StrCallerCudaEngine() {
    free_workspace();
    if (stream_) { cudaStreamDestroy(stream_); stream_ = nullptr; }
}

void StrCallerCudaEngine::allocate_workspace(size_t max_loci, size_t max_reads) {
    free_workspace();
    max_loci_ = max_loci;
    max_reads_ = max_reads;

    cudaMalloc(&d_loci_, max_loci_ * sizeof(StrLocusDescriptor));
    cudaMalloc(&d_reads_, max_reads_ * sizeof(StrReadEvidence));
    cudaMalloc(&d_offsets_, max_loci_ * sizeof(uint32_t));
    cudaMalloc(&d_counts_, max_loci_ * sizeof(uint32_t));
    cudaMalloc(&d_calls_, max_loci_ * sizeof(StrGenotypeCall));
}

void StrCallerCudaEngine::free_workspace() {
    if (d_loci_) { cudaFree(d_loci_); d_loci_ = nullptr; }
    if (d_reads_) { cudaFree(d_reads_); d_reads_ = nullptr; }
    if (d_offsets_) { cudaFree(d_offsets_); d_offsets_ = nullptr; }
    if (d_counts_) { cudaFree(d_counts_); d_counts_ = nullptr; }
    if (d_calls_) { cudaFree(d_calls_); d_calls_ = nullptr; }
}

bool StrCallerCudaEngine::genotype_loci(
    const std::vector<StrLocusDescriptor>& loci,
    const std::vector<StrReadEvidence>& reads,
    const std::vector<uint32_t>& locus_read_offsets,
    const std::vector<uint32_t>& locus_read_counts,
    std::vector<StrGenotypeCall>& out_calls,
    StrExecutionStats& out_stats,
    uint16_t max_repeat_search
) {
    uint32_t num_loci = static_cast<uint32_t>(loci.size());
    size_t num_reads = reads.size();
    if (num_loci == 0) return false;

    if (num_loci > max_loci_ || num_reads > max_reads_) {
        allocate_workspace(num_loci * 2, num_reads * 2);
    }

    auto t_start = std::chrono::high_resolution_clock::now();

    cudaMemcpyAsync(d_loci_, loci.data(), num_loci * sizeof(StrLocusDescriptor), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_reads_, reads.data(), num_reads * sizeof(StrReadEvidence), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_offsets_, locus_read_offsets.data(), num_loci * sizeof(uint32_t), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_counts_, locus_read_counts.data(), num_loci * sizeof(uint32_t), cudaMemcpyHostToDevice, stream_);

    cudaEvent_t evt_start, evt_stop;
    cudaEventCreate(&evt_start);
    cudaEventCreate(&evt_stop);
    cudaEventRecord(evt_start, stream_);

    int threads = 256;
    xoos_str_genotype_solver_kernel<<<num_loci, threads, 0, stream_>>>(
        d_loci_, d_reads_, d_offsets_, d_counts_, num_loci, max_repeat_search, d_calls_
    );

    cudaEventRecord(evt_stop, stream_);

    out_calls.resize(num_loci);
    cudaMemcpyAsync(out_calls.data(), d_calls_, num_loci * sizeof(StrGenotypeCall), cudaMemcpyDeviceToHost, stream_);
    cudaStreamSynchronize(stream_);

    float kernel_time = 0.0f;
    cudaEventElapsedTime(&kernel_time, evt_start, evt_stop);
    cudaEventDestroy(evt_start);
    cudaEventDestroy(evt_stop);

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    out_stats.total_loci_processed = num_loci;
    out_stats.total_reads_evaluated = num_reads;
    out_stats.total_genotypes_tested = static_cast<uint64_t>(num_loci) * ((max_repeat_search * (max_repeat_search + 1)) / 2);
    out_stats.kernel_time_ms = kernel_time;
    out_stats.total_time_ms = total_ms;
    if (kernel_time > 0.0f) {
        out_stats.throughput_reads_per_sec = num_reads / (kernel_time / 1000.0);
    }

    return true;
}

void StrCallerCudaEngine::print_device_info() const {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, device_id_ < count ? device_id_ : 0);
        std::cout << "\n================================================================================" << std::endl;
        std::cout << "  XOOS STR & REPEAT EXPANSION CALLER CUDA ENGINE (sm_" << prop.major << prop.minor << ")" << std::endl;
        std::cout << "================================================================================" << std::endl;
        std::cout << "  Device Name:             " << prop.name << std::endl;
        std::cout << "  Compute Capability:      " << prop.major << "." << prop.minor << std::endl;
        std::cout << "  Streaming Multiprocessors: " << prop.multiProcessorCount << " SMs" << std::endl;
        std::cout << "  Total Global Memory:     " << (prop.totalGlobalMem >> 20) << " MB" << std::endl;
        std::cout << "  L2 Cache Size:           " << (prop.l2CacheSize >> 20) << " MB" << std::endl;
        std::cout << "================================================================================\n" << std::endl;
    }
}

} // namespace xoos::str_caller::cuda
