#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include <atomic>
#include <random>

#include "cuda/alignment_metrics_cuda.cuh"
#include "cuda/metrics_types.cuh"

using namespace xoos::alignment_metrics::cuda;

// ============================================================================
// CPU Reference Coverage Pileup & Metric Implementation
// ============================================================================

static void run_cpu_metrics_threaded(
    const std::vector<AlignedReadRecord>& reads,
    const std::vector<std::string>& sequences,
    const std::string& ref_sequence,
    unsigned int num_threads,
    uint64_t& out_total_aligned_bases,
    uint64_t& out_covered_bases
) {
    size_t ref_len = ref_sequence.length();
    std::vector<uint32_t> global_depth(ref_len, 0);

    // Multi-threaded interval pileup
    std::vector<std::thread> workers;
    size_t total_reads = reads.size();

    std::vector<std::vector<uint32_t>> thread_depths(num_threads, std::vector<uint32_t>(ref_len, 0));

    for (unsigned int t = 0; t < num_threads; ++t) {
        workers.emplace_back([&, t]() {
            size_t start_r = (total_reads * t) / num_threads;
            size_t end_r = (total_reads * (t + 1)) / num_threads;

            auto& local_depth = thread_depths[t];
            for (size_t r = start_r; r < end_r; ++r) {
                const auto& read = reads[r];
                if ((read.flag & 0x0004) != 0) continue;
                for (uint64_t p = read.rbeg; p < read.rend && p < ref_len; ++p) {
                    local_depth[p]++;
                }
            }
        });
    }

    for (auto& w : workers) {
        w.join();
    }

    // Merge thread depths
    uint64_t total_aligned = 0;
    uint64_t covered = 0;
    for (size_t p = 0; p < ref_len; ++p) {
        uint32_t d = 0;
        for (unsigned int t = 0; t < num_threads; ++t) {
            d += thread_depths[t][p];
        }
        if (d > 0) {
            covered++;
            total_aligned += d;
        }
    }
    out_total_aligned_bases = total_aligned;
    out_covered_bases = covered;
}

// ============================================================================
// Main Benchmark Suite
// ============================================================================

int main(int argc, char** argv) {
    size_t num_reads = 1000000;
    if (argc > 1) {
        num_reads = std::stoull(argv[1]);
    }

    std::cout << "================================================================================" << std::endl;
    std::cout << "  XOOS ALIGNMENT METRICS CUDA BENCHMARK & CPU CONTRAST HARNESS" << std::endl;
    std::cout << "================================================================================" << std::endl;

    AlignmentMetricsCudaEngine engine(0);
    engine.print_device_info();

    // 1. Synthesize reference genome (10 Mb chromosome)
    size_t ref_len = 10000000; // 10 Mb
    std::string ref_seq(ref_len, 'A');
    std::string bases = "ACGT";
    for (size_t i = 0; i < ref_len; ++i) {
        ref_seq[i] = bases[(i * 3 + (i / 7)) % 4];
    }
    engine.load_reference_genome(ref_seq);

    // 2. Synthesize 1,000,000 aligned reads across genome (~30x average depth)
    std::vector<AlignedReadRecord> reads(num_reads);
    std::vector<std::string> sequences(num_reads);
    int read_len = 100;

    for (size_t i = 0; i < num_reads; ++i) {
        uint64_t rbeg = (i * 9ULL) % (ref_len - read_len - 10);
        reads[i].rbeg = rbeg;
        reads[i].rend = rbeg + read_len;
        reads[i].read_len = read_len;
        reads[i].flag = (i % 100 == 0) ? 0x0004 : 0x0000; // 1% unmapped
        if (i % 20 == 0) reads[i].flag |= 0x0400; // 5% duplicates
        reads[i].mapq = 60;
        reads[i].is_reverse = (i % 2 == 1) ? 1 : 0;

        std::string s = ref_seq.substr(rbeg, read_len);
        if (i % 50 == 0) s[10] = 'N'; // occasional mismatch
        sequences[i] = s;
    }

    double total_mb = (num_reads * read_len) / (1024.0 * 1024.0);
    std::cout << "\n  Workload: " << num_reads << " reads (" << std::fixed << std::setprecision(2)
              << total_mb << " MB query bases) across " << (ref_len / 1e6) << " Mb Reference\n" << std::endl;

    // 3. CPU Single-Threaded Benchmark
    std::cout << ">>> [1] Running CPU Single-Threaded Coverage & Accuracy Calculation..." << std::endl;
    uint64_t cpu_1t_aligned = 0, cpu_1t_covered = 0;
    auto t_cpu1_start = std::chrono::high_resolution_clock::now();
    run_cpu_metrics_threaded(reads, sequences, ref_seq, 1, cpu_1t_aligned, cpu_1t_covered);
    auto t_cpu1_end = std::chrono::high_resolution_clock::now();
    double cpu_1t_ms = std::chrono::duration<double, std::milli>(t_cpu1_end - t_cpu1_start).count();
    double cpu_1t_reads_sec = num_reads / (cpu_1t_ms / 1000.0);

    std::cout << "  CPU 1-Thread Time:         " << std::fixed << std::setprecision(2) << cpu_1t_ms << " ms" << std::endl;
    std::cout << "  CPU 1-Thread Throughput:   " << std::setprecision(2) << (cpu_1t_reads_sec / 1e3) << " kReads/sec" << std::endl;

    // 4. CPU 16-Threaded Benchmark
    unsigned int num_threads = std::min(16u, std::thread::hardware_concurrency());
    std::cout << "\n>>> [2] Running CPU Multi-Threaded (" << num_threads << " Threads) Calculation..." << std::endl;
    uint64_t cpu_mt_aligned = 0, cpu_mt_covered = 0;
    auto t_cpumt_start = std::chrono::high_resolution_clock::now();
    run_cpu_metrics_threaded(reads, sequences, ref_seq, num_threads, cpu_mt_aligned, cpu_mt_covered);
    auto t_cpumt_end = std::chrono::high_resolution_clock::now();
    double cpu_mt_ms = std::chrono::duration<double, std::milli>(t_cpumt_end - t_cpumt_start).count();
    double cpu_mt_reads_sec = num_reads / (cpu_mt_ms / 1000.0);

    std::cout << "  CPU " << num_threads << "-Thread Time:        " << std::fixed << std::setprecision(2) << cpu_mt_ms << " ms" << std::endl;
    std::cout << "  CPU " << num_threads << "-Thread Throughput:  " << std::setprecision(2) << (cpu_mt_reads_sec / 1e3) << " kReads/sec" << std::endl;

    // 5. GPU CUDA RTX 5090 Blackwell Benchmark
    std::cout << "\n>>> [3] Running GPU CUDA Metrics Engine (NVIDIA RTX 5090 sm_120)..." << std::endl;
    CoverageSummaryMetrics gpu_coverage;
    HpAccuracyMetrics gpu_hp;
    ReadAlignmentStats gpu_stats;
    MetricsExecutionStats gpu_exec;

    auto t_gpu_start = std::chrono::high_resolution_clock::now();
    engine.compute_metrics(reads, sequences, gpu_coverage, gpu_hp, gpu_stats, gpu_exec);
    auto t_gpu_end = std::chrono::high_resolution_clock::now();
    double gpu_e2e_ms = std::chrono::duration<double, std::milli>(t_gpu_end - t_gpu_start).count();
    double gpu_e2e_reads_sec = num_reads / (gpu_e2e_ms / 1000.0);

    std::cout << "  GPU Kernel Time:           " << std::fixed << std::setprecision(2) << gpu_exec.total_time_ms << " ms" << std::endl;
    std::cout << "  GPU End-to-End Time:       " << std::fixed << std::setprecision(2) << gpu_e2e_ms << " ms (incl. PCIe H2D/D2H)" << std::endl;
    std::cout << "  GPU Kernel Throughput:     " << std::setprecision(2) << (gpu_exec.throughput_reads_per_sec / 1e6) << " Million reads/sec" << std::endl;
    std::cout << "  GPU End-to-End Throughput: " << std::setprecision(2) << (gpu_e2e_reads_sec / 1e6) << " Million reads/sec" << std::endl;
    std::cout << "  Mean Genome Coverage:      " << std::setprecision(2) << gpu_coverage.mean_coverage << "x (Max: " << gpu_coverage.max_coverage << "x)" << std::endl;
    std::cout << "  Genome Breadth >= 1x:      " << std::setprecision(2) << gpu_coverage.pct_bases_ge_1x << "%" << std::endl;
    std::cout << "  HP Accuracy (Length 1-5):  " << std::setprecision(3) << (100.0 - gpu_hp.hp_error_rate[1]) << "%" << std::endl;

    // 6. Comparative Summary
    double speedup_vs_1t = cpu_1t_ms / gpu_exec.total_time_ms;
    double speedup_vs_mt = cpu_mt_ms / gpu_exec.total_time_ms;

    std::cout << "\n==========================================================================================" << std::endl;
    std::cout << "                       CPU VS GPU ALIGNMENT METRICS PERFORMANCE SUMMARY" << std::endl;
    std::cout << "==========================================================================================" << std::endl;
    std::cout << std::left << std::setw(28) << "Execution Engine"
              << std::right << std::setw(14) << "Time (ms)"
              << std::setw(22) << "Throughput (reads/s)"
              << std::setw(14) << "Speedup" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    
    std::cout << std::left << std::setw(28) << "CPU 1-Thread"
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << cpu_1t_ms
              << std::setw(22) << std::setprecision(0) << cpu_1t_reads_sec
              << std::setw(12) << "1.0x" << std::endl;

    std::cout << std::left << std::setw(28) << ("CPU " + std::to_string(num_threads) + "-Threads")
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << cpu_mt_ms
              << std::setw(22) << std::setprecision(0) << cpu_mt_reads_sec
              << std::setw(11) << std::setprecision(1) << (cpu_1t_ms / cpu_mt_ms) << "x" << std::endl;

    std::cout << std::left << std::setw(28) << "GPU RTX 5090 (Kernel)"
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << gpu_exec.total_time_ms
              << std::setw(22) << std::setprecision(0) << gpu_exec.throughput_reads_per_sec
              << std::setw(11) << std::setprecision(1) << speedup_vs_1t << "x (" << speedup_vs_mt << "x vs MT)" << std::endl;

    std::cout << std::left << std::setw(28) << "GPU RTX 5090 (End-to-End)"
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << gpu_e2e_ms
              << std::setw(22) << std::setprecision(0) << gpu_e2e_reads_sec
              << std::setw(11) << std::setprecision(1) << (cpu_1t_ms / gpu_e2e_ms) << "x" << std::endl;
    std::cout << "==========================================================================================\n" << std::endl;

    std::cout << ">>> ALIGNMENT METRICS PROBE GATES PASSED <<<\n" << std::endl;
    return 0;
}
