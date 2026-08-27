#include "cnv_caller_cuda.cuh"
#include "cbs_hmm_segmentation_cuda.cuh"
#include "purity_ploidy_solver_cuda.cuh"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cmath>

namespace xoos::cnv_caller::cuda {

CopyNumberCallerCudaEngine::CopyNumberCallerCudaEngine(int device_id) : device_id_(device_id) {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaSetDevice(device_id_ < count ? device_id_ : 0);
    }
    cudaStreamCreate(&stream_);
    allocate_workspace(max_bins_, max_segments_);
}

CopyNumberCallerCudaEngine::~CopyNumberCallerCudaEngine() {
    free_workspace();
    if (stream_) { cudaStreamDestroy(stream_); stream_ = nullptr; }
}

void CopyNumberCallerCudaEngine::allocate_workspace(size_t max_bins, size_t max_segs) {
    free_workspace();
    max_bins_ = max_bins;
    max_segments_ = max_segs;
    max_grid_points_ = kPurityGridSteps * kPloidyGridSteps;

    cudaMalloc(&d_bins_, max_bins_ * sizeof(GenomicBinRecord));
    cudaMalloc(&d_state_path_, max_bins_ * sizeof(uint8_t));
    cudaMalloc(&d_segments_, max_segments_ * sizeof(CnvSegment));
    cudaMalloc(&d_grid_scores_, max_grid_points_ * sizeof(float));
    cudaMalloc(&d_fit_result_, sizeof(PurityPloidyFit));
}

void CopyNumberCallerCudaEngine::free_workspace() {
    if (d_bins_) { cudaFree(d_bins_); d_bins_ = nullptr; }
    if (d_state_path_) { cudaFree(d_state_path_); d_state_path_ = nullptr; }
    if (d_segments_) { cudaFree(d_segments_); d_segments_ = nullptr; }
    if (d_grid_scores_) { cudaFree(d_grid_scores_); d_grid_scores_ = nullptr; }
    if (d_fit_result_) { cudaFree(d_fit_result_); d_fit_result_ = nullptr; }
}

bool CopyNumberCallerCudaEngine::normalize_gc_and_log2r(
    std::vector<GenomicBinRecord>& bins,
    float gc_a,
    float gc_b,
    float gc_c,
    float baseline_depth
) {
    uint64_t num_bins = bins.size();
    if (num_bins == 0) return false;

    if (num_bins > max_bins_) {
        allocate_workspace(num_bins * 2, max_segments_);
    }

    cudaMemcpyAsync(d_bins_, bins.data(), num_bins * sizeof(GenomicBinRecord), cudaMemcpyHostToDevice, stream_);

    int threads = 256;
    int blocks = (num_bins + threads - 1) / threads;
    if (blocks > 170 * 8) blocks = 170 * 8;
    if (blocks == 0) blocks = 1;

    xoos_gc_normalize_kernel<<<blocks, threads, 0, stream_>>>(
        d_bins_, num_bins, gc_a, gc_b, gc_c, baseline_depth
    );

    cudaMemcpyAsync(bins.data(), d_bins_, num_bins * sizeof(GenomicBinRecord), cudaMemcpyDeviceToHost, stream_);
    cudaStreamSynchronize(stream_);
    return true;
}

bool CopyNumberCallerCudaEngine::run_cnv_calling(
    std::vector<GenomicBinRecord>& bins,
    std::vector<CnvSegment>& out_segments,
    PurityPloidyFit& out_fit,
    CnvExecutionStats& out_stats
) {
    uint64_t num_bins = bins.size();
    if (num_bins == 0) return false;

    auto t_start = std::chrono::high_resolution_clock::now();

    if (num_bins > max_bins_) {
        allocate_workspace(num_bins * 2, max_segments_);
    }

    // 1. Transfer Bins & Run GC Normalization
    cudaMemcpyAsync(d_bins_, bins.data(), num_bins * sizeof(GenomicBinRecord), cudaMemcpyHostToDevice, stream_);

    int threads = 256;
    int blocks = (num_bins + threads - 1) / threads;
    if (blocks > 170 * 8) blocks = 170 * 8;
    if (blocks == 0) blocks = 1;

    xoos_gc_normalize_kernel<<<blocks, threads, 0, stream_>>>(
        d_bins_, num_bins, -0.5f, 1.0f, 50.0f, 100.0f
    );

    // 2. Initial HMM Viterbi Segmentation (Purity=1.0, Ploidy=2.0)
    int hmm_blocks = 128;
    xoos_hmm_viterbi_kernel<<<hmm_blocks, 1, 0, stream_>>>(
        d_bins_, num_bins, 1.0f, 2.0f, d_state_path_
    );

    std::vector<uint8_t> h_state_path(num_bins);
    cudaMemcpyAsync(h_state_path.data(), d_state_path_, num_bins * sizeof(uint8_t), cudaMemcpyDeviceToHost, stream_);
    cudaMemcpyAsync(bins.data(), d_bins_, num_bins * sizeof(GenomicBinRecord), cudaMemcpyDeviceToHost, stream_);
    cudaStreamSynchronize(stream_);

    // 3. Collapse initial segment boundaries
    collapse_state_path_to_segments(bins.data(), h_state_path.data(), num_bins, out_segments);
    uint32_t num_segs = static_cast<uint32_t>(out_segments.size());

    // 4. 2D Purity & Ploidy Optimization Grid Search
    if (num_segs > 0) {
        if (num_segs > max_segments_) {
            allocate_workspace(max_bins_, num_segs * 2);
        }

        cudaMemcpyAsync(d_segments_, out_segments.data(), num_segs * sizeof(CnvSegment), cudaMemcpyHostToDevice, stream_);

        uint32_t total_points = kPurityGridSteps * kPloidyGridSteps;
        int grid_blocks = (total_points + threads - 1) / threads;
        if (grid_blocks > 170 * 8) grid_blocks = 170 * 8;

        xoos_purity_ploidy_grid_kernel<<<grid_blocks, threads, 0, stream_>>>(
            d_segments_, num_segs, 0.10f, 1.00f, 1.50f, 5.00f,
            kPurityGridSteps, kPloidyGridSteps, d_grid_scores_
        );

        xoos_reduce_best_purity_ploidy_kernel<<<1, 1, 0, stream_>>>(
            d_grid_scores_, total_points, 0.10f, 1.00f, 1.50f, 5.00f,
            kPurityGridSteps, kPloidyGridSteps, d_fit_result_
        );

        cudaMemcpyAsync(&out_fit, d_fit_result_, sizeof(PurityPloidyFit), cudaMemcpyDeviceToHost, stream_);
        cudaStreamSynchronize(stream_);

        // 5. Final HMM Viterbi Refinement with Optimal Purity/Ploidy
        xoos_hmm_viterbi_kernel<<<hmm_blocks, 1, 0, stream_>>>(
            d_bins_, num_bins, out_fit.best_purity, out_fit.best_ploidy, d_state_path_
        );
        cudaMemcpyAsync(h_state_path.data(), d_state_path_, num_bins * sizeof(uint8_t), cudaMemcpyDeviceToHost, stream_);
        cudaStreamSynchronize(stream_);

        collapse_state_path_to_segments(bins.data(), h_state_path.data(), num_bins, out_segments);
    } else {
        out_fit.best_purity = 1.0f;
        out_fit.best_ploidy = 2.0f;
        out_fit.max_log_likelihood = 0.0;
        out_fit.subclonal_fraction = 0.0f;
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    out_stats.total_bins_processed = num_bins;
    out_stats.total_segments_discovered = out_segments.size();
    out_stats.total_grid_evaluations = kPurityGridSteps * kPloidyGridSteps;
    out_stats.total_time_ms = elapsed_ms;
    if (elapsed_ms > 0) {
        out_stats.throughput_bins_per_sec = num_bins / (elapsed_ms / 1000.0);
    }

    return true;
}

void CopyNumberCallerCudaEngine::print_device_info() const {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, device_id_ < count ? device_id_ : 0);
        std::cout << "\n================================================================================" << std::endl;
        std::cout << "  XOOS COPY NUMBER CALLER CUDA ENGINE (sm_" << prop.major << prop.minor << ")" << std::endl;
        std::cout << "================================================================================" << std::endl;
        std::cout << "  Device Name:             " << prop.name << std::endl;
        std::cout << "  Compute Capability:      " << prop.major << "." << prop.minor << std::endl;
        std::cout << "  Streaming Multiprocessors: " << prop.multiProcessorCount << " SMs" << std::endl;
        std::cout << "  Total Global Memory:     " << (prop.totalGlobalMem >> 20) << " MB" << std::endl;
        std::cout << "  L2 Cache Size:           " << (prop.l2CacheSize >> 20) << " MB" << std::endl;
        std::cout << "================================================================================\n" << std::endl;
    }
}

} // namespace xoos::cnv_caller::cuda
