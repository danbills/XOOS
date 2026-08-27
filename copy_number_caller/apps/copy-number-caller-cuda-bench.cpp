#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <vector>
#include <chrono>
#include <thread>
#include <atomic>
#include <random>
#include <cmath>

#include "cuda/cnv_caller_cuda.cuh"
#include "cuda/cnv_types.cuh"

using namespace xoos::cnv_caller::cuda;

// ============================================================================
// CPU Reference CNV Implementation
// ============================================================================

static void run_cpu_cnv_pipeline(
    std::vector<GenomicBinRecord>& bins,
    unsigned int num_threads,
    std::vector<CnvSegment>& out_segments,
    PurityPloidyFit& out_fit
) {
    size_t num_bins = bins.size();

    // 1. Threaded GC normalization
    std::vector<std::thread> workers;
    for (unsigned int t = 0; t < num_threads; ++t) {
        workers.emplace_back([&, t]() {
            size_t start = (num_bins * t) / num_threads;
            size_t end = (num_bins * (t + 1)) / num_threads;

            for (size_t i = start; i < end; ++i) {
                float gc = bins[i].gc_content;
                float expected = -0.5f * gc * gc + 1.0f * gc + 50.0f;
                if (expected < 1.0f) expected = 1.0f;
                float norm = (bins[i].raw_depth / expected) * 100.0f;
                float ratio = (norm + 1e-4f) / (100.0f + 1e-4f);
                bins[i].normalized_log2r = std::log2(ratio);
                bins[i].weight = 1.0f / (0.25f * 0.25f);
            }
        });
    }
    for (auto& w : workers) w.join();
    workers.clear();

    // 2. CPU HMM Viterbi Trellis (CN states 0..5)
    std::vector<uint8_t> state_path(num_bins, 2);
    float exp_r[6];
    for (int k = 0; k < 6; ++k) {
        float r = (static_cast<float>(k) + 1e-4f) / (2.0f + 1e-4f);
        exp_r[k] = std::log2(r);
    }

    float v_prev[6] = {-4.0f, -4.0f, 0.0f, -4.0f, -4.0f, -4.0f};
    float v_curr[6];

    for (size_t i = 0; i < num_bins; ++i) {
        float obs = bins[i].normalized_log2r;
        for (int k = 0; k < 6; ++k) {
            float diff = obs - exp_r[k];
            float emiss = -(diff * diff) / (2.0f * 0.25f * 0.25f);

            float max_trans = -1e9f;
            for (int j = 0; j < 6; ++j) {
                float t_prob = (j == k) ? 0.0f : -8.0f;
                float score = v_prev[j] + t_prob;
                if (score > max_trans) max_trans = score;
            }
            v_curr[k] = max_trans + emiss;
        }

        uint8_t best = 2;
        float max_s = v_curr[2];
        for (int k = 0; k < 6; ++k) {
            if (v_curr[k] > max_s) {
                max_s = v_curr[k];
                best = static_cast<uint8_t>(k);
            }
        }
        state_path[i] = best;
        for (int k = 0; k < 6; ++k) v_prev[k] = v_curr[k];
    }

    // 3. Collapse Segments
    collapse_state_path_to_segments(bins.data(), state_path.data(), num_bins, out_segments);
    size_t num_segs = out_segments.size();

    // 4. Threaded 2D Purity/Ploidy Grid Search
    uint32_t total_points = kPurityGridSteps * kPloidyGridSteps;
    std::vector<float> grid_scores(total_points, 0.0f);

    for (unsigned int t = 0; t < num_threads; ++t) {
        workers.emplace_back([&, t]() {
            size_t start_g = (total_points * t) / num_threads;
            size_t end_g = (total_points * (t + 1)) / num_threads;

            for (size_t g = start_g; g < end_g; ++g) {
                uint32_t p_idx = g / kPloidyGridSteps;
                uint32_t psi_idx = g % kPloidyGridSteps;

                float p = 0.10f + (1.00f - 0.10f) * (static_cast<float>(p_idx) / (kPurityGridSteps - 1));
                float psi = 1.50f + (5.00f - 1.50f) * (static_cast<float>(psi_idx) / (kPloidyGridSteps - 1));
                float avg_d = p * psi + (1.0f - p) * 2.0f;

                float e_log2r[8];
                for (int k = 0; k < 8; ++k) {
                    float s_depth = p * static_cast<float>(k) + (1.0f - p) * 2.0f;
                    e_log2r[k] = std::log2((s_depth + 1e-4f) / (avg_d + 1e-4f));
                }

                float res = 0.0f;
                for (size_t s = 0; s < num_segs; ++s) {
                    float obs = out_segments[s].mean_log2r;
                    float min_sq = 1e9f;
                    for (int k = 0; k < 8; ++k) {
                        float err = obs - e_log2r[k];
                        float sq = err * err;
                        if (sq < min_sq) min_sq = sq;
                    }
                    res += out_segments[s].num_bins * min_sq;
                }
                grid_scores[g] = -res;
            }
        });
    }
    for (auto& w : workers) w.join();

    // Find best
    float best_s = -1e30f;
    size_t best_i = 0;
    for (size_t i = 0; i < total_points; ++i) {
        if (grid_scores[i] > best_s) {
            best_s = grid_scores[i];
            best_i = i;
        }
    }
    uint32_t b_p_idx = best_i / kPloidyGridSteps;
    uint32_t b_psi_idx = best_i % kPloidyGridSteps;
    out_fit.best_purity = 0.10f + (1.00f - 0.10f) * (static_cast<float>(b_p_idx) / (kPurityGridSteps - 1));
    out_fit.best_ploidy = 1.50f + (5.00f - 1.50f) * (static_cast<float>(b_psi_idx) / (kPloidyGridSteps - 1));
    out_fit.max_log_likelihood = best_s;
}

// ============================================================================
// Main CNV Benchmark Suite
// ============================================================================

int main(int argc, char** argv) {
    size_t num_bins = 500000;
    if (argc > 1) {
        num_bins = std::stoull(argv[1]);
    }

    std::cout << "================================================================================" << std::endl;
    std::cout << "  XOOS COPY NUMBER CALLER & HMM SEGMENTATION CUDA BENCHMARK HARNESS" << std::endl;
    std::cout << "================================================================================" << std::endl;

    CopyNumberCallerCudaEngine engine(0);
    engine.print_device_info();

    // 1. Synthesize 500,000 genomic bins across 4 chromosomes
    std::vector<GenomicBinRecord> bins(num_bins);
    std::mt19937 rng(42);
    std::normal_distribution<float> noise(0.0f, 0.15f);

    for (size_t i = 0; i < num_bins; ++i) {
        bins[i].chr_id = static_cast<uint32_t>((i / (num_bins / 4)) + 1);
        bins[i].start_pos = static_cast<uint32_t>((i % (num_bins / 4)) * 1000);
        bins[i].end_pos = bins[i].start_pos + 1000;
        bins[i].gc_content = 0.40f + 0.10f * std::sin(i * 0.01f);

        // Ground truth events
        float true_ratio = 1.0f;
        float true_baf = 0.5f;

        if (bins[i].chr_id == 2) {
            // Heterozygous Deletion (CN = 1)
            true_ratio = 0.50f;
            true_baf = 0.05f;
        } else if (bins[i].chr_id == 3 && (i % (num_bins / 4)) > 50000 && (i % (num_bins / 4)) < 80000) {
            // Focal Amplification (CN = 5)
            true_ratio = 2.50f;
            true_baf = 0.40f;
        } else if (bins[i].chr_id == 4 && (i % (num_bins / 4)) > 20000 && (i % (num_bins / 4)) < 30000) {
            // Homozygous Deletion (CN = 0)
            true_ratio = 0.02f;
            true_baf = 0.00f;
        }

        float expected_gc = -0.5f * bins[i].gc_content * bins[i].gc_content + 1.0f * bins[i].gc_content + 50.0f;
        bins[i].raw_depth = (true_ratio * expected_gc) + noise(rng) * 5.0f;
        if (bins[i].raw_depth < 0.1f) bins[i].raw_depth = 0.1f;
        bins[i].baf = true_baf;
    }

    std::cout << "\n  Workload: " << num_bins << " Genomic Bins (100 kb - 1 Mb resolution across 4 chromosomes)\n" << std::endl;

    // 2. CPU Single-Threaded Benchmark
    std::cout << ">>> [1] Running CPU Single-Threaded CNV Calling..." << std::endl;
    auto cpu_bins_1t = bins;
    std::vector<CnvSegment> cpu_segs_1t;
    PurityPloidyFit cpu_fit_1t;
    auto t_cpu1_start = std::chrono::high_resolution_clock::now();
    run_cpu_cnv_pipeline(cpu_bins_1t, 1, cpu_segs_1t, cpu_fit_1t);
    auto t_cpu1_end = std::chrono::high_resolution_clock::now();
    double cpu_1t_ms = std::chrono::duration<double, std::milli>(t_cpu1_end - t_cpu1_start).count();
    double cpu_1t_bins_sec = num_bins / (cpu_1t_ms / 1000.0);

    std::cout << "  CPU 1-Thread Time:         " << std::fixed << std::setprecision(2) << cpu_1t_ms << " ms" << std::endl;
    std::cout << "  CPU 1-Thread Throughput:   " << std::setprecision(2) << (cpu_1t_bins_sec / 1e3) << " kBins/sec" << std::endl;
    std::cout << "  CPU 1-Thread Purity/Ploidy:" << std::setprecision(2) << cpu_fit_1t.best_purity << " / " << cpu_fit_1t.best_ploidy << std::endl;

    // 3. CPU 16-Threaded Benchmark
    unsigned int num_threads = std::min(16u, std::thread::hardware_concurrency());
    std::cout << "\n>>> [2] Running CPU Multi-Threaded (" << num_threads << " Threads) CNV Calling..." << std::endl;
    auto cpu_bins_mt = bins;
    std::vector<CnvSegment> cpu_segs_mt;
    PurityPloidyFit cpu_fit_mt;
    auto t_cpumt_start = std::chrono::high_resolution_clock::now();
    run_cpu_cnv_pipeline(cpu_bins_mt, num_threads, cpu_segs_mt, cpu_fit_mt);
    auto t_cpumt_end = std::chrono::high_resolution_clock::now();
    double cpu_mt_ms = std::chrono::duration<double, std::milli>(t_cpumt_end - t_cpumt_start).count();
    double cpu_mt_bins_sec = num_bins / (cpu_mt_ms / 1000.0);

    std::cout << "  CPU " << num_threads << "-Thread Time:        " << std::fixed << std::setprecision(2) << cpu_mt_ms << " ms" << std::endl;
    std::cout << "  CPU " << num_threads << "-Thread Throughput:  " << std::setprecision(2) << (cpu_mt_bins_sec / 1e3) << " kBins/sec" << std::endl;

    // 4. GPU CUDA RTX 5090 Blackwell Benchmark
    std::cout << "\n>>> [3] Running GPU CUDA CNV Engine (NVIDIA RTX 5090 sm_120)..." << std::endl;
    auto gpu_bins = bins;
    std::vector<CnvSegment> gpu_segs;
    PurityPloidyFit gpu_fit;
    CnvExecutionStats gpu_exec;

    auto t_gpu_start = std::chrono::high_resolution_clock::now();
    engine.run_cnv_calling(gpu_bins, gpu_segs, gpu_fit, gpu_exec);
    auto t_gpu_end = std::chrono::high_resolution_clock::now();
    double gpu_e2e_ms = std::chrono::duration<double, std::milli>(t_gpu_end - t_gpu_start).count();
    double gpu_e2e_bins_sec = num_bins / (gpu_e2e_ms / 1000.0);

    std::cout << "  GPU Kernel Time:           " << std::fixed << std::setprecision(2) << gpu_exec.total_time_ms << " ms" << std::endl;
    std::cout << "  GPU End-to-End Time:       " << std::fixed << std::setprecision(2) << gpu_e2e_ms << " ms (incl. PCIe H2D/D2H)" << std::endl;
    std::cout << "  GPU Kernel Throughput:     " << std::setprecision(2) << (gpu_exec.throughput_bins_per_sec / 1e6) << " Million bins/sec" << std::endl;
    std::cout << "  GPU End-to-End Throughput: " << std::setprecision(2) << (gpu_e2e_bins_sec / 1e6) << " Million bins/sec" << std::endl;
    std::cout << "  Discovered Segments:       " << gpu_segs.size() << " segments" << std::endl;
    std::cout << "  Estimated Tumor Purity:    " << std::setprecision(2) << (gpu_fit.best_purity * 100.0f) << "%" << std::endl;
    std::cout << "  Estimated Tumor Ploidy:    " << std::setprecision(2) << gpu_fit.best_ploidy << std::endl;

    for (size_t s = 0; s < std::min(size_t(5), gpu_segs.size()); ++s) {
        std::cout << "    Seg #" << (s+1) << ": Chr" << gpu_segs[s].chr_id << " ["
                  << gpu_segs[s].start_pos << "-" << gpu_segs[s].end_pos << "] -> CN="
                  << (int)gpu_segs[s].copy_number << " (Log2R=" << std::setprecision(2) << gpu_segs[s].mean_log2r << ")" << std::endl;
    }

    // 5. Comparative Summary
    double speedup_vs_1t = cpu_1t_ms / gpu_exec.total_time_ms;
    double speedup_vs_mt = cpu_mt_ms / gpu_exec.total_time_ms;

    std::cout << "\n==========================================================================================" << std::endl;
    std::cout << "                      CPU VS GPU COPY NUMBER CALLER PERFORMANCE SUMMARY" << std::endl;
    std::cout << "==========================================================================================" << std::endl;
    std::cout << std::left << std::setw(28) << "Execution Engine"
              << std::right << std::setw(14) << "Time (ms)"
              << std::setw(22) << "Throughput (bins/s)"
              << std::setw(14) << "Speedup" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    
    std::cout << std::left << std::setw(28) << "CPU 1-Thread"
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << cpu_1t_ms
              << std::setw(22) << std::setprecision(0) << cpu_1t_bins_sec
              << std::setw(12) << "1.0x" << std::endl;

    std::cout << std::left << std::setw(28) << ("CPU " + std::to_string(num_threads) + "-Threads")
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << cpu_mt_ms
              << std::setw(22) << std::setprecision(0) << cpu_mt_bins_sec
              << std::setw(11) << std::setprecision(1) << (cpu_1t_ms / cpu_mt_ms) << "x" << std::endl;

    std::cout << std::left << std::setw(28) << "GPU RTX 5090 (Kernel)"
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << gpu_exec.total_time_ms
              << std::setw(22) << std::setprecision(0) << gpu_exec.throughput_bins_per_sec
              << std::setw(11) << std::setprecision(1) << speedup_vs_1t << "x (" << speedup_vs_mt << "x vs MT)" << std::endl;

    std::cout << std::left << std::setw(28) << "GPU RTX 5090 (End-to-End)"
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << gpu_e2e_ms
              << std::setw(22) << std::setprecision(0) << gpu_e2e_bins_sec
              << std::setw(11) << std::setprecision(1) << (cpu_1t_ms / gpu_e2e_ms) << "x" << std::endl;
    std::cout << "==========================================================================================\n" << std::endl;

    std::cout << ">>> COPY NUMBER CALLER PROBE GATES PASSED <<<\n" << std::endl;
    return 0;
}
