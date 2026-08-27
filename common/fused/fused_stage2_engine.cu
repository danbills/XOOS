#include "fused_stage2_engine.cuh"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cmath>

namespace xoos::fused::cuda {

FusedStage2CudaEngine::FusedStage2CudaEngine(int device_id, size_t initial_max_reads) 
    : device_id_(device_id), max_reads_(initial_max_reads) {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaSetDevice(device_id_ < count ? device_id_ : 0);
    }
    cudaStreamCreate(&stream_);
    allocate_workspace(max_reads_);
    jit_compiler_ = std::make_unique<jit::NvrtcSuperKernelJit>();
}

FusedStage2CudaEngine::~FusedStage2CudaEngine() {
    free_workspace();
    if (stream_) { cudaStreamDestroy(stream_); stream_ = nullptr; }
}

void FusedStage2CudaEngine::allocate_workspace(size_t max_reads) {
    free_workspace();
    max_reads_ = max_reads;

    cudaMalloc(&d_reads_, max_reads_ * sizeof(FusedReadRecord));
    cudaMallocHost(&h_pinned_reads_, max_reads_ * sizeof(FusedReadRecord));
    cudaMalloc(&d_metrics_, sizeof(GlobalMetricsAccumulator));
}

void FusedStage2CudaEngine::free_workspace() {
    if (d_reads_) { cudaFree(d_reads_); d_reads_ = nullptr; }
    if (h_pinned_reads_) { cudaFreeHost(h_pinned_reads_); h_pinned_reads_ = nullptr; }
    if (d_metrics_) { cudaFree(d_metrics_); d_metrics_ = nullptr; }
}

FusedReadRecord* FusedStage2CudaEngine::get_pinned_host_buffer(size_t required_reads) {
    if (required_reads > max_reads_) {
        allocate_workspace(required_reads * 2);
    }
    return h_pinned_reads_;
}

bool FusedStage2CudaEngine::execute_aot_canonical(
    std::vector<FusedReadRecord>& reads,
    GlobalMetricsAccumulator& out_metrics,
    FusedExecutionStats& out_stats
) {
    uint64_t num_reads = reads.size();
    if (num_reads == 0) return false;

    if (num_reads > max_reads_) {
        allocate_workspace(num_reads * 2);
    }

    auto t_start = std::chrono::high_resolution_clock::now();

    cudaMemsetAsync(d_metrics_, 0, sizeof(GlobalMetricsAccumulator), stream_);
    cudaMemcpyAsync(d_reads_, reads.data(), num_reads * sizeof(FusedReadRecord), cudaMemcpyHostToDevice, stream_);

    cudaEvent_t evt_start, evt_stop;
    cudaEventCreate(&evt_start);
    cudaEventCreate(&evt_stop);
    cudaEventRecord(evt_start, stream_);

    int threads = 256;
    int blocks = (num_reads + threads - 1) / threads;
    if (blocks > 170 * 8) blocks = 170 * 8;
    if (blocks == 0) blocks = 1;

    xoos_fused_stage2_super_kernel<Canonical_DeepCfDna_Policy><<<blocks, threads, 0, stream_>>>(
        d_reads_, num_reads, d_metrics_
    );

    cudaEventRecord(evt_stop, stream_);

    cudaMemcpyAsync(&out_metrics, d_metrics_, sizeof(GlobalMetricsAccumulator), cudaMemcpyDeviceToHost, stream_);
    cudaMemcpyAsync(reads.data(), d_reads_, num_reads * sizeof(FusedReadRecord), cudaMemcpyDeviceToHost, stream_);
    cudaStreamSynchronize(stream_);

    float kernel_time = 0.0f;
    cudaEventElapsedTime(&kernel_time, evt_start, evt_stop);
    cudaEventDestroy(evt_start);
    cudaEventDestroy(evt_stop);

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    out_stats.total_reads_processed = num_reads;
    out_stats.jit_compile_time_ms = 0.0;
    out_stats.kernel_time_ms = kernel_time;
    out_stats.total_time_ms = total_ms;
    out_stats.used_jit = false;
    if (kernel_time > 0.0f) {
        out_stats.throughput_reads_per_sec = num_reads / (kernel_time / 1000.0);
        double bytes_processed = static_cast<double>(num_reads * sizeof(FusedReadRecord) * 2);
        out_stats.vram_bandwidth_gb_per_sec = (bytes_processed / 1e9) / (kernel_time / 1000.0);
    }

    return true;
}

bool FusedStage2CudaEngine::execute_jit_dynamic(
    const jit::DynamicPolicyConfig& config,
    std::vector<FusedReadRecord>& reads,
    GlobalMetricsAccumulator& out_metrics,
    FusedExecutionStats& out_stats
) {
    uint64_t num_reads = reads.size();
    if (num_reads == 0) return false;

    if (num_reads > max_reads_) {
        allocate_workspace(num_reads * 2);
    }

    CUmodule module = nullptr;
    CUfunction func = nullptr;
    double compile_ms = 0.0;

    bool jit_ok = jit_compiler_->get_or_compile_kernel(config, module, func, compile_ms);
    if (!jit_ok || !func) return false;

    auto t_start = std::chrono::high_resolution_clock::now();

    cudaMemsetAsync(d_metrics_, 0, sizeof(GlobalMetricsAccumulator), stream_);
    cudaMemcpyAsync(d_reads_, reads.data(), num_reads * sizeof(FusedReadRecord), cudaMemcpyHostToDevice, stream_);

    cudaEvent_t evt_start, evt_stop;
    cudaEventCreate(&evt_start);
    cudaEventCreate(&evt_stop);
    cudaEventRecord(evt_start, stream_);

    int threads = 256;
    int blocks = (num_reads + threads - 1) / threads;
    if (blocks > 170 * 8) blocks = 170 * 8;
    if (blocks == 0) blocks = 1;

    void* args[] = { &d_reads_, &num_reads, &d_metrics_ };

    bool launch_ok = jit_compiler_->launch_kernel(
        func,
        blocks,
        threads,
        reinterpret_cast<CUstream>(stream_),
        args
    );

    if (!launch_ok) return false;

    cudaEventRecord(evt_stop, stream_);

    cudaMemcpyAsync(&out_metrics, d_metrics_, sizeof(GlobalMetricsAccumulator), cudaMemcpyDeviceToHost, stream_);
    cudaMemcpyAsync(reads.data(), d_reads_, num_reads * sizeof(FusedReadRecord), cudaMemcpyDeviceToHost, stream_);
    cudaStreamSynchronize(stream_);

    float kernel_time = 0.0f;
    cudaEventElapsedTime(&kernel_time, evt_start, evt_stop);
    cudaEventDestroy(evt_start);
    cudaEventDestroy(evt_stop);

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    out_stats.total_reads_processed = num_reads;
    out_stats.jit_compile_time_ms = compile_ms;
    out_stats.kernel_time_ms = kernel_time;
    out_stats.total_time_ms = total_ms;
    out_stats.used_jit = true;
    if (kernel_time > 0.0f) {
        out_stats.throughput_reads_per_sec = num_reads / (kernel_time / 1000.0);
        double bytes_processed = static_cast<double>(num_reads * sizeof(FusedReadRecord) * 2);
        out_stats.vram_bandwidth_gb_per_sec = (bytes_processed / 1e9) / (kernel_time / 1000.0);
    }

    return true;
}

bool FusedStage2CudaEngine::execute_jit_dynamic_pinned(
    const jit::DynamicPolicyConfig& config,
    FusedReadRecord* h_pinned_reads,
    uint64_t num_reads,
    GlobalMetricsAccumulator& out_metrics,
    FusedExecutionStats& out_stats
) {
    if (!h_pinned_reads || num_reads == 0) return false;

    if (num_reads > max_reads_) {
        allocate_workspace(num_reads * 2);
    }

    CUmodule module = nullptr;
    CUfunction func = nullptr;
    double compile_ms = 0.0;

    bool jit_ok = jit_compiler_->get_or_compile_kernel(config, module, func, compile_ms);
    if (!jit_ok || !func) return false;

    auto t_start = std::chrono::high_resolution_clock::now();

    cudaMemsetAsync(d_metrics_, 0, sizeof(GlobalMetricsAccumulator), stream_);
    // Direct zero-copy PCIe DMA transfer
    cudaMemcpyAsync(d_reads_, h_pinned_reads, num_reads * sizeof(FusedReadRecord), cudaMemcpyHostToDevice, stream_);

    cudaEvent_t evt_start, evt_stop;
    cudaEventCreate(&evt_start);
    cudaEventCreate(&evt_stop);
    cudaEventRecord(evt_start, stream_);

    int threads = 256;
    int blocks = (num_reads + threads - 1) / threads;
    if (blocks > 170 * 8) blocks = 170 * 8;
    if (blocks == 0) blocks = 1;

    void* args[] = { &d_reads_, &num_reads, &d_metrics_ };

    bool launch_ok = jit_compiler_->launch_kernel(
        func,
        blocks,
        threads,
        reinterpret_cast<CUstream>(stream_),
        args
    );

    if (!launch_ok) return false;

    cudaEventRecord(evt_stop, stream_);

    cudaMemcpyAsync(&out_metrics, d_metrics_, sizeof(GlobalMetricsAccumulator), cudaMemcpyDeviceToHost, stream_);
    cudaMemcpyAsync(h_pinned_reads, d_reads_, num_reads * sizeof(FusedReadRecord), cudaMemcpyDeviceToHost, stream_);
    cudaStreamSynchronize(stream_);

    float kernel_time = 0.0f;
    cudaEventElapsedTime(&kernel_time, evt_start, evt_stop);
    cudaEventDestroy(evt_start);
    cudaEventDestroy(evt_stop);

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    out_stats.total_reads_processed = num_reads;
    out_stats.jit_compile_time_ms = compile_ms;
    out_stats.kernel_time_ms = kernel_time;
    out_stats.total_time_ms = total_ms;
    out_stats.used_jit = true;
    if (kernel_time > 0.0f) {
        out_stats.throughput_reads_per_sec = num_reads / (kernel_time / 1000.0);
        double bytes_processed = static_cast<double>(num_reads * sizeof(FusedReadRecord) * 2);
        out_stats.vram_bandwidth_gb_per_sec = (bytes_processed / 1e9) / (kernel_time / 1000.0);
    }

    return true;
}

void FusedStage2CudaEngine::print_device_info() const {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, device_id_ < count ? device_id_ : 0);
        std::cout << "\n================================================================================" << std::endl;
        std::cout << "  XOOS FUSED STAGE-2 SUPER-KERNEL ENGINE (sm_" << prop.major << prop.minor << ")" << std::endl;
        std::cout << "================================================================================" << std::endl;
        std::cout << "  Device Name:             " << prop.name << std::endl;
        std::cout << "  Compute Capability:      " << prop.major << "." << prop.minor << std::endl;
        std::cout << "  Streaming Multiprocessors: " << prop.multiProcessorCount << " SMs" << std::endl;
        std::cout << "  Total Global Memory:     " << (prop.totalGlobalMem >> 20) << " MB" << std::endl;
        std::cout << "  L2 Cache Size:           " << (prop.l2CacheSize >> 20) << " MB" << std::endl;
        std::cout << "================================================================================\n" << std::endl;
    }
}

} // namespace xoos::fused::cuda
