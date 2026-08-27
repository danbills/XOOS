#include "tumor_fraction_cuda.cuh"
#include "variant_pileup_filter_cuda.cuh"
#include "noise_correction_solver_cuda.cuh"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cmath>

namespace xoos::tfe::cuda {

TumorFractionCudaEngine::TumorFractionCudaEngine(int device_id) : device_id_(device_id) {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaSetDevice(device_id_ < count ? device_id_ : 0);
    }
    cudaStreamCreate(&stream_);
    allocate_workspace(max_sites_, max_reads_);
}

TumorFractionCudaEngine::~TumorFractionCudaEngine() {
    free_workspace();
    if (stream_) { cudaStreamDestroy(stream_); stream_ = nullptr; }
}

void TumorFractionCudaEngine::allocate_workspace(size_t max_sites, size_t max_reads) {
    free_workspace();
    max_sites_ = max_sites;
    max_reads_ = max_reads;

    cudaMalloc(&d_sites_, max_sites_ * sizeof(VariantProbeSite));
    cudaMalloc(&d_reads_, max_reads_ * sizeof(ProbeReadObservation));
    cudaMalloc(&d_offsets_, max_sites_ * sizeof(uint32_t));
    cudaMalloc(&d_counts_, max_sites_ * sizeof(uint32_t));
    cudaMalloc(&d_results_, max_sites_ * sizeof(VariantSitePileupResult));
    cudaMalloc(&d_summary_, sizeof(TumorFractionSummary));
}

void TumorFractionCudaEngine::free_workspace() {
    if (d_sites_) { cudaFree(d_sites_); d_sites_ = nullptr; }
    if (d_reads_) { cudaFree(d_reads_); d_reads_ = nullptr; }
    if (d_offsets_) { cudaFree(d_offsets_); d_offsets_ = nullptr; }
    if (d_counts_) { cudaFree(d_counts_); d_counts_ = nullptr; }
    if (d_results_) { cudaFree(d_results_); d_results_ = nullptr; }
    if (d_summary_) { cudaFree(d_summary_); d_summary_ = nullptr; }
}

bool TumorFractionCudaEngine::estimate_tumor_fraction(
    const std::vector<VariantProbeSite>& sites,
    const std::vector<ProbeReadObservation>& reads,
    const std::vector<uint32_t>& site_offsets,
    const std::vector<uint32_t>& site_counts,
    std::vector<VariantSitePileupResult>& out_results,
    TumorFractionSummary& out_summary,
    TfeExecutionStats& out_stats
) {
    uint32_t num_sites = static_cast<uint32_t>(sites.size());
    size_t num_reads = reads.size();
    if (num_sites == 0) return false;

    if (num_sites > max_sites_ || num_reads > max_reads_) {
        allocate_workspace(num_sites * 2, num_reads * 2);
    }

    auto t_start = std::chrono::high_resolution_clock::now();

    cudaMemcpyAsync(d_sites_, sites.data(), num_sites * sizeof(VariantProbeSite), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_reads_, reads.data(), num_reads * sizeof(ProbeReadObservation), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_offsets_, site_offsets.data(), num_sites * sizeof(uint32_t), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_counts_, site_counts.data(), num_sites * sizeof(uint32_t), cudaMemcpyHostToDevice, stream_);

    cudaEvent_t evt_start, evt_stop;
    cudaEventCreate(&evt_start);
    cudaEventCreate(&evt_stop);
    cudaEventRecord(evt_start, stream_);

    int threads = 256;
    xoos_variant_pileup_kernel<<<num_sites, threads, 0, stream_>>>(
        d_sites_, d_reads_, d_offsets_, d_counts_, num_sites, d_results_
    );

    xoos_tfe_summary_solver_kernel<<<1, 1, 0, stream_>>>(
        d_results_, num_sites, d_summary_
    );

    cudaEventRecord(evt_stop, stream_);

    out_results.resize(num_sites);
    cudaMemcpyAsync(out_results.data(), d_results_, num_sites * sizeof(VariantSitePileupResult), cudaMemcpyDeviceToHost, stream_);
    cudaMemcpyAsync(&out_summary, d_summary_, sizeof(TumorFractionSummary), cudaMemcpyDeviceToHost, stream_);
    cudaStreamSynchronize(stream_);

    float kernel_time = 0.0f;
    cudaEventElapsedTime(&kernel_time, evt_start, evt_stop);
    cudaEventDestroy(evt_start);
    cudaEventDestroy(evt_stop);

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    out_stats.total_sites_evaluated = num_sites;
    out_stats.total_reads_processed = num_reads;
    out_stats.pileup_time_ms = kernel_time;
    out_stats.solver_time_ms = kernel_time;
    out_stats.total_time_ms = total_ms;
    if (kernel_time > 0.0f) {
        out_stats.throughput_reads_per_sec = num_reads / (kernel_time / 1000.0);
    }

    return true;
}

void TumorFractionCudaEngine::print_device_info() const {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, device_id_ < count ? device_id_ : 0);
        std::cout << "\n================================================================================" << std::endl;
        std::cout << "  XOOS TUMOR FRACTION ESTIMATOR CUDA ENGINE (sm_" << prop.major << prop.minor << ")" << std::endl;
        std::cout << "================================================================================" << std::endl;
        std::cout << "  Device Name:             " << prop.name << std::endl;
        std::cout << "  Compute Capability:      " << prop.major << "." << prop.minor << std::endl;
        std::cout << "  Streaming Multiprocessors: " << prop.multiProcessorCount << " SMs" << std::endl;
        std::cout << "  Total Global Memory:     " << (prop.totalGlobalMem >> 20) << " MB" << std::endl;
        std::cout << "  L2 Cache Size:           " << (prop.l2CacheSize >> 20) << " MB" << std::endl;
        std::cout << "================================================================================\n" << std::endl;
    }
}

} // namespace xoos::tfe::cuda
