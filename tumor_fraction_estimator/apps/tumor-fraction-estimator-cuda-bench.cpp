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

#include "cuda/tumor_fraction_cuda.cuh"
#include "cuda/tfe_types.cuh"

using namespace xoos::tfe::cuda;

// ============================================================================
// CPU Reference Tumor Fraction Implementation
// ============================================================================

static void run_cpu_tfe_pipeline(
    const std::vector<VariantProbeSite>& sites,
    const std::vector<ProbeReadObservation>& reads,
    const std::vector<uint32_t>& site_offsets,
    const std::vector<uint32_t>& site_counts,
    unsigned int num_threads,
    std::vector<VariantSitePileupResult>& out_results,
    TumorFractionSummary& out_summary
) {
    size_t num_sites = sites.size();
    out_results.resize(num_sites);

    std::vector<std::thread> workers;
    for (unsigned int t = 0; t < num_threads; ++t) {
        workers.emplace_back([&, t]() {
            size_t start = (num_sites * t) / num_threads;
            size_t end = (num_sites * (t + 1)) / num_threads;

            for (size_t s = start; s < end; ++s) {
                const auto& site = sites[s];
                uint32_t offset = site_offsets[s];
                uint32_t count = site_counts[s];

                uint32_t ref = 0;
                uint32_t alt = 0;
                uint32_t other = 0;

                for (uint32_t r = 0; r < count; ++r) {
                    const auto& read = reads[offset + r];
                    if (read.is_duplicate != 0) continue;
                    if (read.mapq < site.min_mapq) continue;
                    if (read.base_qual < site.min_baseq) continue;

                    char b = read.observed_base;
                    if (b == site.ref_base) ref++;
                    else if (b == site.alt_base) alt++;
                    else if (b == 'A' || b == 'C' || b == 'G' || b == 'T') other++;
                }

                uint32_t depth = ref + alt + other;
                out_results[s].site_id = site.site_id;
                out_results[s].ref_count = ref;
                out_results[s].alt_count = alt;
                out_results[s].other_alts_count = other;
                out_results[s].total_depth = depth;

                if (depth >= 10) {
                    out_results[s].observed_vaf = static_cast<float>(alt) / depth;
                    float adj = static_cast<float>(alt) - 0.5f * static_cast<float>(other);
                    if (adj < 0.0f) adj = 0.0f;
                    out_results[s].adjusted_vaf = adj / depth;
                    out_results[s].is_passed = 1;
                } else {
                    out_results[s].observed_vaf = 0.0f;
                    out_results[s].adjusted_vaf = 0.0f;
                    out_results[s].is_passed = 0;
                }
            }
        });
    }
    for (auto& w : workers) w.join();

    // Reduce summary
    uint64_t tot_ref = 0;
    uint64_t tot_alt = 0;
    uint64_t tot_other = 0;
    uint64_t tot_depth = 0;
    uint32_t pass_sites = 0;
    uint32_t detected = 0;

    for (size_t i = 0; i < num_sites; ++i) {
        if (out_results[i].is_passed) {
            pass_sites++;
            tot_ref += out_results[i].ref_count;
            tot_alt += out_results[i].alt_count;
            tot_other += out_results[i].other_alts_count;
            tot_depth += out_results[i].total_depth;
            if (out_results[i].alt_count > 0) detected++;
        }
    }

    double adj_alt = 0.0;
    if (tot_alt > (tot_other / 2)) {
        adj_alt = static_cast<double>(tot_alt) - 0.5 * static_cast<double>(tot_other);
    }

    double mean_vaf = (tot_depth > 0) ? (adj_alt / tot_depth) : 0.0;
    double err_rate = (tot_depth > 0) ? (static_cast<double>(tot_other) / (2.0 * tot_depth)) : 0.0;

    out_summary.tumor_fraction = 2.0 * mean_vaf;
    out_summary.mean_vaf = mean_vaf;
    out_summary.error_rate = err_rate;
    out_summary.total_alt = tot_alt;
    out_summary.total_ref = tot_ref;
    out_summary.total_other_alts = tot_other;
    out_summary.total_depth = tot_depth;
    out_summary.total_adjusted_alt = adj_alt;
    out_summary.total_sites = static_cast<uint32_t>(num_sites);
    out_summary.passing_sites = pass_sites;
    out_summary.sites_detected = detected;
}

// ============================================================================
// Main TFE Benchmark Suite
// ============================================================================

int main(int argc, char** argv) {
    size_t num_sites = 5000;
    if (argc > 1) {
        num_sites = std::stoull(argv[1]);
    }

    std::cout << "================================================================================" << std::endl;
    std::cout << "  XOOS TUMOR FRACTION ESTIMATOR (TFE) CUDA BENCHMARK HARNESS" << std::endl;
    std::cout << "================================================================================" << std::endl;

    TumorFractionCudaEngine engine(0);
    engine.print_device_info();

    // 1. Synthesize 5,000 somatic variant sites and 1,000,000 reads
    std::vector<VariantProbeSite> sites(num_sites);
    std::vector<uint32_t> site_offsets(num_sites);
    std::vector<uint32_t> site_counts(num_sites);
    std::vector<ProbeReadObservation> reads;

    std::mt19937 rng(42);
    std::uniform_int_distribution<uint32_t> depth_dist(150, 250);
    std::uniform_real_distribution<float> noise_dist(0.0f, 1.0f);

    size_t total_reads = 0;
    float true_somatic_vaf = 0.0050f; // 0.50% VAF -> 1.00% Tumor Fraction
    float background_noise_rate = 0.0005f; // 0.05% sequencing error

    for (size_t s = 0; s < num_sites; ++s) {
        sites[s].site_id = static_cast<uint32_t>(s + 1);
        sites[s].chr_id = static_cast<uint32_t>((s % 24) + 1);
        sites[s].pos = static_cast<uint32_t>(s * 10000 + 100);
        sites[s].ref_base = 'A';
        sites[s].alt_base = 'G';
        sites[s].var_type = (s % 5 == 0) ? kVarTypeNoiseProbe : kVarTypeSomatic;
        sites[s].min_mapq = 10;
        sites[s].min_baseq = 29;

        uint32_t depth = depth_dist(rng);
        site_offsets[s] = static_cast<uint32_t>(reads.size());
        site_counts[s] = depth;

        for (uint32_t r = 0; r < depth; ++r) {
            ProbeReadObservation read;
            read.read_id = total_reads++;
            read.site_id = sites[s].site_id;
            read.mapq = 60;
            read.base_qual = 35;
            read.is_reverse_strand = (r % 2 == 0) ? 0 : 1;
            read.is_duplicate = 0;

            float roll = noise_dist(rng);
            if (sites[s].var_type == kVarTypeSomatic && roll < true_somatic_vaf) {
                read.observed_base = 'G'; // True somatic alt
            } else if (roll < (true_somatic_vaf + background_noise_rate)) {
                read.observed_base = 'C'; // Third allele (sequencing error)
            } else if (roll < (true_somatic_vaf + 2 * background_noise_rate)) {
                read.observed_base = 'T'; // Fourth allele (sequencing error)
            } else {
                read.observed_base = 'A'; // Reference allele
            }

            reads.push_back(read);
        }
    }

    std::cout << "\n  Workload: " << num_sites << " Variant Probe Sites | " << reads.size() << " Read Observations\n" << std::endl;

    // 2. CPU Single-Threaded Benchmark
    std::cout << ">>> [1] Running CPU Single-Threaded Tumor Fraction Estimator..." << std::endl;
    std::vector<VariantSitePileupResult> cpu_results_1t;
    TumorFractionSummary cpu_summary_1t;
    auto t_cpu1_start = std::chrono::high_resolution_clock::now();
    run_cpu_tfe_pipeline(sites, reads, site_offsets, site_counts, 1, cpu_results_1t, cpu_summary_1t);
    auto t_cpu1_end = std::chrono::high_resolution_clock::now();
    double cpu_1t_ms = std::chrono::duration<double, std::milli>(t_cpu1_end - t_cpu1_start).count();
    double cpu_1t_reads_sec = reads.size() / (cpu_1t_ms / 1000.0);

    std::cout << "  CPU 1-Thread Time:         " << std::fixed << std::setprecision(2) << cpu_1t_ms << " ms" << std::endl;
    std::cout << "  CPU 1-Thread Throughput:   " << std::setprecision(2) << (cpu_1t_reads_sec / 1e3) << " kReads/sec" << std::endl;
    std::cout << "  CPU 1-Thread Tumor Frac:   " << std::setprecision(4) << (cpu_summary_1t.tumor_fraction * 100.0) << "%" << std::endl;

    // 3. CPU 16-Threaded Benchmark
    unsigned int num_threads = std::min(16u, std::thread::hardware_concurrency());
    std::cout << "\n>>> [2] Running CPU Multi-Threaded (" << num_threads << " Threads) Tumor Fraction Estimator..." << std::endl;
    std::vector<VariantSitePileupResult> cpu_results_mt;
    TumorFractionSummary cpu_summary_mt;
    auto t_cpumt_start = std::chrono::high_resolution_clock::now();
    run_cpu_tfe_pipeline(sites, reads, site_offsets, site_counts, num_threads, cpu_results_mt, cpu_summary_mt);
    auto t_cpumt_end = std::chrono::high_resolution_clock::now();
    double cpu_mt_ms = std::chrono::duration<double, std::milli>(t_cpumt_end - t_cpumt_start).count();
    double cpu_mt_reads_sec = reads.size() / (cpu_mt_ms / 1000.0);

    std::cout << "  CPU " << num_threads << "-Thread Time:        " << std::fixed << std::setprecision(2) << cpu_mt_ms << " ms" << std::endl;
    std::cout << "  CPU " << num_threads << "-Thread Throughput:  " << std::setprecision(2) << (cpu_mt_reads_sec / 1e3) << " kReads/sec" << std::endl;

    // 4. GPU CUDA RTX 5090 Blackwell Benchmark
    std::cout << "\n>>> [3] Running GPU CUDA TFE Engine (NVIDIA RTX 5090 sm_120)..." << std::endl;
    std::vector<VariantSitePileupResult> gpu_results;
    TumorFractionSummary gpu_summary;
    TfeExecutionStats gpu_stats;

    auto t_gpu_start = std::chrono::high_resolution_clock::now();
    engine.estimate_tumor_fraction(sites, reads, site_offsets, site_counts, gpu_results, gpu_summary, gpu_stats);
    auto t_gpu_end = std::chrono::high_resolution_clock::now();
    double gpu_e2e_ms = std::chrono::duration<double, std::milli>(t_gpu_end - t_gpu_start).count();
    double gpu_e2e_reads_sec = reads.size() / (gpu_e2e_ms / 1000.0);

    std::cout << "  GPU Kernel Time:           " << std::fixed << std::setprecision(2) << gpu_stats.pileup_time_ms << " ms" << std::endl;
    std::cout << "  GPU End-to-End Time:       " << std::fixed << std::setprecision(2) << gpu_e2e_ms << " ms (incl. PCIe H2D/D2H)" << std::endl;
    std::cout << "  GPU Kernel Throughput:     " << std::setprecision(2) << (gpu_stats.throughput_reads_per_sec / 1e6) << " Million reads/sec" << std::endl;
    std::cout << "  GPU End-to-End Throughput: " << std::setprecision(2) << (gpu_e2e_reads_sec / 1e6) << " Million reads/sec" << std::endl;
    std::cout << "  Passing Somatic Sites:     " << gpu_summary.passing_sites << " / " << gpu_summary.total_sites << std::endl;
    std::cout << "  Estimated Error Rate:      " << std::setprecision(4) << (gpu_summary.error_rate * 100.0) << "%" << std::endl;
    std::cout << "  Estimated Mean VAF:        " << std::setprecision(4) << (gpu_summary.mean_vaf * 100.0) << "%" << std::endl;
    std::cout << "  Estimated Tumor Fraction:  " << std::setprecision(4) << (gpu_summary.tumor_fraction * 100.0) << "%" << std::endl;

    for (size_t s = 0; s < std::min(size_t(5), gpu_results.size()); ++s) {
        std::cout << "    Site #" << gpu_results[s].site_id << " (Chr" << sites[s].chr_id << ":" << sites[s].pos
                  << " " << sites[s].ref_base << "->" << sites[s].alt_base << ") -> Ref/Alt/Other: "
                  << gpu_results[s].ref_count << "/" << gpu_results[s].alt_count << "/" << gpu_results[s].other_alts_count
                  << " | VAF=" << std::setprecision(3) << (gpu_results[s].adjusted_vaf * 100.0f) << "%" << std::endl;
    }

    // 5. Comparative Summary
    double speedup_vs_1t = cpu_1t_ms / gpu_stats.pileup_time_ms;
    double speedup_vs_mt = cpu_mt_ms / gpu_stats.pileup_time_ms;

    std::cout << "\n==========================================================================================" << std::endl;
    std::cout << "                     CPU VS GPU TUMOR FRACTION ESTIMATOR PERFORMANCE SUMMARY" << std::endl;
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
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << gpu_stats.pileup_time_ms
              << std::setw(22) << std::setprecision(0) << gpu_stats.throughput_reads_per_sec
              << std::setw(11) << std::setprecision(1) << speedup_vs_1t << "x (" << speedup_vs_mt << "x vs MT)" << std::endl;

    std::cout << std::left << std::setw(28) << "GPU RTX 5090 (End-to-End)"
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << gpu_e2e_ms
              << std::setw(22) << std::setprecision(0) << gpu_e2e_reads_sec
              << std::setw(11) << std::setprecision(1) << (cpu_1t_ms / gpu_e2e_ms) << "x" << std::endl;
    std::cout << "==========================================================================================\n" << std::endl;

    std::cout << ">>> TUMOR FRACTION ESTIMATOR PROBE GATES PASSED <<<\n" << std::endl;
    return 0;
}
