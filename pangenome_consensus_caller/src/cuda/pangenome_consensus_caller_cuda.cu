#include "pangenome_consensus_caller_cuda.cuh"
#include "yc_decoder_cuda.cuh"
#include "pangenome_rescuer_cuda.cuh"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cmath>

namespace xoos::pangenome::cuda {

PangenomeConsensusCallerCudaEngine::PangenomeConsensusCallerCudaEngine(int device_id) : device_id_(device_id) {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaSetDevice(device_id_ < count ? device_id_ : 0);
    }
    cudaStreamCreate(&stream_);
    allocate_workspace(max_reads_);
}

PangenomeConsensusCallerCudaEngine::~PangenomeConsensusCallerCudaEngine() {
    free_workspace();
    if (stream_) { cudaStreamDestroy(stream_); stream_ = nullptr; }
}

void PangenomeConsensusCallerCudaEngine::allocate_workspace(size_t max_reads) {
    free_workspace();
    max_reads_ = max_reads;

    cudaMalloc(&d_reads_, max_reads_ * sizeof(DuplexReadRecord));
    cudaMalloc(&d_results_, max_reads_ * sizeof(PangenomeUpdateResult));
}

void PangenomeConsensusCallerCudaEngine::free_workspace() {
    if (d_reads_) { cudaFree(d_reads_); d_reads_ = nullptr; }
    if (d_results_) { cudaFree(d_results_); d_results_ = nullptr; }
}

bool PangenomeConsensusCallerCudaEngine::rescue_consensus_reads(
    std::vector<DuplexReadRecord>& reads,
    std::vector<PangenomeUpdateResult>& out_results,
    PangenomeExecutionStats& out_stats
) {
    uint64_t num_reads = reads.size();
    if (num_reads == 0) return false;

    if (num_reads > max_reads_) {
        allocate_workspace(num_reads * 2);
    }

    auto t_start = std::chrono::high_resolution_clock::now();

    cudaMemcpyAsync(d_reads_, reads.data(), num_reads * sizeof(DuplexReadRecord), cudaMemcpyHostToDevice, stream_);

    cudaEvent_t evt_start, evt_stop;
    cudaEventCreate(&evt_start);
    cudaEventCreate(&evt_stop);
    cudaEventRecord(evt_start, stream_);

    int threads = 256;
    int blocks = (num_reads + threads - 1) / threads;
    if (blocks > 170 * 8) blocks = 170 * 8;
    if (blocks == 0) blocks = 1;

    xoos_pangenome_rescue_kernel<<<blocks, threads, 0, stream_>>>(
        d_reads_, num_reads, d_results_
    );

    cudaEventRecord(evt_stop, stream_);

    out_results.resize(num_reads);
    cudaMemcpyAsync(out_results.data(), d_results_, num_reads * sizeof(PangenomeUpdateResult), cudaMemcpyDeviceToHost, stream_);
    cudaMemcpyAsync(reads.data(), d_reads_, num_reads * sizeof(DuplexReadRecord), cudaMemcpyDeviceToHost, stream_);
    cudaStreamSynchronize(stream_);

    float kernel_time = 0.0f;
    cudaEventElapsedTime(&kernel_time, evt_start, evt_stop);
    cudaEventDestroy(evt_start);
    cudaEventDestroy(evt_stop);

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    uint64_t corrections = 0;
    for (const auto& r : out_results) {
        corrections += r.num_corrections;
    }

    out_stats.total_reads_processed = num_reads;
    out_stats.total_corrections_made = corrections;
    out_stats.total_discrepancies_evaluated = num_reads * 150;
    out_stats.kernel_time_ms = kernel_time;
    out_stats.total_time_ms = total_ms;
    if (kernel_time > 0.0f) {
        out_stats.throughput_reads_per_sec = num_reads / (kernel_time / 1000.0);
    }

    return true;
}

void PangenomeConsensusCallerCudaEngine::print_device_info() const {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, device_id_ < count ? device_id_ : 0);
        std::cout << "\n================================================================================" << std::endl;
        std::cout << "  XOOS PANGENOME CONSENSUS CALLER CUDA ENGINE (sm_" << prop.major << prop.minor << ")" << std::endl;
        std::cout << "================================================================================" << std::endl;
        std::cout << "  Device Name:             " << prop.name << std::endl;
        std::cout << "  Compute Capability:      " << prop.major << "." << prop.minor << std::endl;
        std::cout << "  Streaming Multiprocessors: " << prop.multiProcessorCount << " SMs" << std::endl;
        std::cout << "  Total Global Memory:     " << (prop.totalGlobalMem >> 20) << " MB" << std::endl;
        std::cout << "  L2 Cache Size:           " << (prop.l2CacheSize >> 20) << " MB" << std::endl;
        std::cout << "================================================================================\n" << std::endl;
    }
}

} // namespace xoos::pangenome::cuda
