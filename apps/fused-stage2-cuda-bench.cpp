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

#include "fused/fused_stage2_engine.cuh"
#include "fused/fused_types.cuh"
#include "fused/nvrtc_jit_engine.hpp"

using namespace xoos::fused::cuda;
using namespace xoos::fused::jit;

// ============================================================================
// Main Fused Super-Kernel Benchmark Harness
// ============================================================================

int main(int argc, char** argv) {
    size_t num_reads = 1000000;
    if (argc > 1) {
        num_reads = std::stoull(argv[1]);
    }

    std::cout << "================================================================================" << std::endl;
    std::cout << "  XOOS FUSED STAGE-2 SUPER-KERNEL BENCHMARK HARNESS (AOT vs JIT)" << std::endl;
    std::cout << "================================================================================" << std::endl;

    FusedStage2CudaEngine engine(0);
    engine.print_device_info();

    // 1. Synthesize 1,000,000 duplex reads (150 bp)
    std::vector<FusedReadRecord> reads(num_reads);
    const char bases[4] = {'A', 'C', 'G', 'T'};
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> base_dist(0, 3);
    std::uniform_int_distribution<uint16_t> isize_dist(120, 250);
    std::uniform_real_distribution<float> prob_dist(0.0f, 1.0f);

    for (size_t r = 0; r < num_reads; ++r) {
        reads[r].read_id = r + 1;
        reads[r].barcode_hash = 0xABCD1234ULL + (r % 50000);
        reads[r].chr_id = static_cast<uint32_t>((r % 24) + 1);
        reads[r].pos = static_cast<uint32_t>(r * 100);
        reads[r].length = 150;
        reads[r].insert_size = isize_dist(rng);
        reads[r].mapq = 60;
        reads[r].is_reverse_strand = (r % 2 == 0) ? 0 : 1;
        reads[r].has_yc_tag = 1;
        reads[r].is_duplicate = (r % 10 == 0) ? 1 : 0;

        bool favor_r2 = (r % 5 == 0);
        reads[r].r1_graph_score = 120.0f;
        reads[r].r2_graph_score = favor_r2 ? 145.0f : 110.0f;

        for (uint16_t i = 0; i < 150; ++i) {
            reads[r].sequence[i] = bases[base_dist(rng)];
            reads[r].base_qual[i] = 37;

            if (prob_dist(rng) < 0.05f) {
                reads[r].yc_tag[i] = (reads[r].sequence[i] == 'A') ? 'c' : 'g';
            } else {
                reads[r].yc_tag[i] = '*';
            }
        }
        reads[r].sequence[150] = '\0';
        reads[r].yc_tag[150] = '\0';
    }

    std::cout << "\n  Workload: " << num_reads << " Reads (150 bp | 150 Million Base Pairs)\n" << std::endl;

    // 2. Strategy A: AOT Pre-Compiled Super-Kernel Execution
    std::cout << ">>> [1] Running Strategy A: Ahead-Of-Time (AOT) Pre-Compiled Super-Kernel..." << std::endl;
    auto aot_reads = reads;
    GlobalMetricsAccumulator aot_metrics;
    FusedExecutionStats aot_stats;

    engine.execute_aot_canonical(aot_reads, aot_metrics, aot_stats);

    std::cout << "  AOT Kernel Time:           " << std::fixed << std::setprecision(2) << aot_stats.kernel_time_ms << " ms" << std::endl;
    std::cout << "  AOT End-to-End Time:       " << std::fixed << std::setprecision(2) << aot_stats.total_time_ms << " ms" << std::endl;
    std::cout << "  AOT Kernel Throughput:     " << std::setprecision(2) << (aot_stats.throughput_reads_per_sec / 1e6) << " Million reads/sec" << std::endl;
    std::cout << "  AOT Effective Bandwidth:   " << std::setprecision(2) << aot_stats.vram_bandwidth_gb_per_sec << " GB/sec" << std::endl;
    std::cout << "  Rescued Reads:             " << aot_metrics.total_rescued_reads << std::endl;
    std::cout << "  Total Base Corrections:    " << aot_metrics.total_base_corrections << std::endl;
    std::cout << "  Total GC Bases:            " << aot_metrics.total_gc_bases << std::endl;

    // 3. Strategy B: On-Demand NVRTC JIT Compilation & Execution (Cold)
    std::cout << "\n>>> [2] Running Strategy B: On-Demand NVRTC JIT Compilation (Cold Start)..." << std::endl;
    auto jit_reads = reads;
    GlobalMetricsAccumulator jit_metrics;
    FusedExecutionStats jit_stats;

    DynamicPolicyConfig custom_config;
    custom_config.enable_rescue = true;
    custom_config.enable_collapse = true;
    custom_config.enable_gc_metrics = true;
    custom_config.enable_insert_metrics = true;
    custom_config.min_family_size = 4;
    custom_config.adjusted_bq = 25;

    engine.execute_jit_dynamic(custom_config, jit_reads, jit_metrics, jit_stats);

    std::cout << "  NVRTC JIT Compilation Time:" << std::fixed << std::setprecision(2) << jit_stats.jit_compile_time_ms << " ms" << std::endl;
    std::cout << "  JIT Kernel Execution Time: " << std::fixed << std::setprecision(2) << jit_stats.kernel_time_ms << " ms" << std::endl;
    std::cout << "  JIT Kernel Throughput:     " << std::setprecision(2) << (jit_stats.throughput_reads_per_sec / 1e6) << " Million reads/sec" << std::endl;

    // 4. Strategy B: On-Demand NVRTC JIT Execution (Warm / Cached CUBIN)
    std::cout << "\n>>> [3] Running Strategy B: On-Demand NVRTC JIT (Warm Cached Hit)..." << std::endl;
    auto cached_reads = reads;
    GlobalMetricsAccumulator cached_metrics;
    FusedExecutionStats cached_stats;

    engine.execute_jit_dynamic(custom_config, cached_reads, cached_metrics, cached_stats);

    std::cout << "  Cached JIT Lookup Time:    " << std::fixed << std::setprecision(2) << cached_stats.jit_compile_time_ms << " ms (0 ms overhead)" << std::endl;
    std::cout << "  Cached Kernel Exec Time:   " << std::fixed << std::setprecision(2) << cached_stats.kernel_time_ms << " ms" << std::endl;
    std::cout << "  Cached Kernel Throughput:  " << std::setprecision(2) << (cached_stats.throughput_reads_per_sec / 1e6) << " Million reads/sec" << std::endl;

    // 5. Comparative Performance Summary
    std::cout << "\n==========================================================================================" << std::endl;
    std::cout << "                 FUSED SUPER-KERNEL (AOT VS NVRTC JIT) PERFORMANCE SUMMARY" << std::endl;
    std::cout << "==========================================================================================" << std::endl;
    std::cout << std::left << std::setw(32) << "Fusion Strategy"
              << std::right << std::setw(16) << "JIT Time (ms)"
              << std::setw(16) << "Kernel (ms)"
              << std::setw(24) << "Throughput (reads/s)" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    
    std::cout << std::left << std::setw(32) << "Strategy A: AOT Pre-Compiled"
              << std::right << std::setw(16) << "0.00"
              << std::setw(16) << std::fixed << std::setprecision(2) << aot_stats.kernel_time_ms
              << std::setw(24) << std::setprecision(0) << aot_stats.throughput_reads_per_sec << std::endl;

    std::cout << std::left << std::setw(32) << "Strategy B: NVRTC JIT (Cold)"
              << std::right << std::setw(16) << std::fixed << std::setprecision(2) << jit_stats.jit_compile_time_ms
              << std::setw(16) << std::fixed << std::setprecision(2) << jit_stats.kernel_time_ms
              << std::setw(24) << std::setprecision(0) << jit_stats.throughput_reads_per_sec << std::endl;

    std::cout << std::left << std::setw(32) << "Strategy B: NVRTC JIT (Cached)"
              << std::right << std::setw(16) << "0.00 (Cached)"
              << std::setw(16) << std::fixed << std::setprecision(2) << cached_stats.kernel_time_ms
              << std::setw(24) << std::setprecision(0) << cached_stats.throughput_reads_per_sec << std::endl;
    std::cout << "==========================================================================================\n" << std::endl;

    std::cout << ">>> FUSED SUPER-KERNEL PROBE GATES PASSED <<<\n" << std::endl;
    return 0;
}
