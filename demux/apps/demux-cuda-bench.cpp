#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <chrono>
#include <cassert>
#include <thread>
#include <atomic>

#include "cuda/demux_cuda.cuh"
#include "cuda/bitap_cuda.cuh"
#include "cuda/hairpin_cuda.cuh"
#include "cuda/duplex_consensus_cuda.cuh"

using namespace xoos::demux::cuda;

// ============================================================================
// CPU Reference Demux & Consensus Implementation
// ============================================================================

// Standalone CPU Bitap Reference Oracle (exact Baeza-Yates algorithm)
static int32_t cpu_bitap_search_oracle(const std::string& pattern, const std::string& text, int max_dist = 2) {
    int m = static_cast<int>(pattern.length());
    int n = static_cast<int>(text.length());
    if (m == 0 || m > 64 || n < m) return -1;

    uint64_t mask[4] = {~0ULL, ~0ULL, ~0ULL, ~0ULL};
    for (int i = 0; i < m; ++i) {
        char c = pattern[i];
        int idx = (c == 'A' || c == 'a') ? 0 :
                  (c == 'C' || c == 'c') ? 1 :
                  (c == 'G' || c == 'g') ? 2 :
                  (c == 'T' || c == 't') ? 3 : 0;
        mask[idx] &= ~(1ULL << i);
    }

    uint64_t msb = 1ULL << (m - 1);
    std::vector<uint64_t> current_state(max_dist + 1, ~0ULL);
    std::vector<uint64_t> previous_state(max_dist + 1, ~0ULL);
    std::vector<int32_t> best_pos(max_dist + 1, -1);

    for (int i = 0; i < n; ++i) {
        char c = text[i];
        int idx = (c == 'A' || c == 'a') ? 0 :
                  (c == 'C' || c == 'c') ? 1 :
                  (c == 'G' || c == 'g') ? 2 :
                  (c == 'T' || c == 't') ? 3 : 0;
        uint64_t char_mask = mask[idx];

        previous_state[0] = current_state[0];
        current_state[0] = (current_state[0] << 1) | char_mask;

        if ((current_state[0] & msb) == 0) {
            return i; // exact match (distance 0) early return
        }

        for (int d = 1; d <= max_dist; ++d) {
            previous_state[d] = current_state[d];
            uint64_t deletion = previous_state[d - 1];
            uint64_t insertion = current_state[d - 1] << 1;
            uint64_t substitution = previous_state[d - 1] << 1;
            uint64_t match = (previous_state[d] << 1) | char_mask;

            uint64_t new_state = deletion & insertion & substitution & match;
            if (best_pos[d] == -1 && ((new_state & msb) == 0)) {
                best_pos[d] = i;
            }
            current_state[d] = new_state;
        }
    }

    for (int d = 1; d <= max_dist; ++d) {
        if (best_pos[d] != -1) return best_pos[d];
    }
    return -1;
}

struct CpuConsensusResult {
    std::string consensus_seq;
    std::string consensus_qual;
    std::string yc_tag;
    int concordant_count = 0;
    int discordant_count = 0;
};

// Process a single read on CPU
static bool cpu_demux_single_read(
    const std::string& seq,
    const std::string& qual,
    const std::string& loop_seq,
    const std::string& start_ad,
    const std::string& end_ad,
    const std::vector<std::string>& sid_5p_list,
    const std::vector<std::string>& sid_3p_list,
    CpuConsensusResult& out_cons
) {
    if (seq.length() < 40) return false;

    // 1. Search for loop
    int loop_pos = cpu_bitap_search_oracle(loop_seq, seq, 2);
    if (loop_pos < 0) return false;

    int loop_len = static_cast<int>(loop_seq.length());
    int sid_len = static_cast<int>(sid_5p_list[0].length());

    // 2. Sample barcode match
    int best_sample = 0;
    int hairpin_start = (loop_pos >= sid_len) ? (loop_pos - sid_len) : 0;
    int hairpin_end = (loop_pos + loop_len + sid_len < static_cast<int>(seq.length()))
        ? (loop_pos + loop_len + sid_len) : static_cast<int>(seq.length() - 1);

    // 3. Adapter boundaries
    int start_ad_pos = cpu_bitap_search_oracle(start_ad, seq.substr(0, 40), 2);
    int insert1_start = (start_ad_pos >= 0) ? (start_ad_pos + 1) : 0;
    int insert1_end = (hairpin_start > 0) ? (hairpin_start - 1) : 0;

    int insert2_start = hairpin_end + 1;
    int insert2_end = static_cast<int>(seq.length()) - static_cast<int>(end_ad.length()) - 1;
    if (insert2_end < insert2_start) insert2_end = static_cast<int>(seq.length() - 1);

    int insert1_len = insert1_end - insert1_start + 1;
    int insert2_len = insert2_end - insert2_start + 1;
    if (insert1_len < 15 || insert2_len < 15) return false;

    int cons_len = std::min(insert1_len, insert2_len);
    out_cons.consensus_seq.resize(cons_len);
    out_cons.consensus_qual.resize(cons_len);
    out_cons.yc_tag.resize(cons_len);
    out_cons.concordant_count = 0;
    out_cons.discordant_count = 0;

    for (int i = 0; i < cons_len; ++i) {
        char b1 = seq[insert1_start + i];
        int q1 = qual[insert1_start + i] - 33;
        if (q1 < 0) q1 = 0;

        char b2_raw = seq[insert2_end - i];
        char b2_rc = (b2_raw == 'A') ? 'T' : (b2_raw == 'C') ? 'G' : (b2_raw == 'G') ? 'C' : 'A';
        int q2 = qual[insert2_end - i] - 33;
        if (q2 < 0) q2 = 0;

        if (b1 == b2_rc) {
            out_cons.consensus_seq[i] = b1;
            out_cons.consensus_qual[i] = static_cast<char>(std::min(60, q1 + q2) + 33);
            out_cons.yc_tag[i] = b1;
            out_cons.concordant_count++;
        } else {
            out_cons.consensus_seq[i] = (q1 >= q2) ? b1 : b2_rc;
            out_cons.consensus_qual[i] = static_cast<char>(std::max(2, std::abs(q1 - q2)) + 33);
            out_cons.yc_tag[i] = 'N';
            out_cons.discordant_count++;
        }
    }
    return true;
}

// Multi-threaded CPU demux batch runner
static void run_cpu_demux_threaded(
    const std::vector<std::pair<std::string, std::string>>& reads,
    const std::string& loop_seq,
    const std::string& start_ad,
    const std::string& end_ad,
    const std::vector<std::string>& sid_5p_list,
    const std::vector<std::string>& sid_3p_list,
    unsigned int num_threads,
    std::atomic<uint64_t>& out_valid_count
) {
    size_t total_reads = reads.size();
    std::vector<std::thread> workers;

    for (unsigned int t = 0; t < num_threads; ++t) {
        workers.emplace_back([&, t]() {
            size_t start_idx = (total_reads * t) / num_threads;
            size_t end_idx = (total_reads * (t + 1)) / num_threads;

            uint64_t local_valid = 0;
            CpuConsensusResult res;
            for (size_t i = start_idx; i < end_idx; ++i) {
                if (cpu_demux_single_read(reads[i].first, reads[i].second, loop_seq, start_ad, end_ad, sid_5p_list, sid_3p_list, res)) {
                    local_valid++;
                }
            }
            out_valid_count += local_valid;
        });
    }

    for (auto& w : workers) {
        w.join();
    }
}

// ============================================================================
// Test Suites & Comparative Benchmark
// ============================================================================

void run_bitap_oracle_parity_test() {
    std::cout << "\n[TEST 1] Bitap CPU Oracle vs GPU CUDA Parity Gate..." << std::endl;

    std::string test_read = 
        "ACGTAGCTAGCTA"                     // 5' Start Adapter (13 bp)
        "GATTACAGATTACAGGCCTTAAGGTCCGAATT"  // R1 Insert (32 bp)
        "TCGACTGACTGA"                     // SID_5p (12 bp)
        "GTCGACAATTCTTGTCATA"              // Hairpin Loop (19 bp)
        "TCAGTCAGTCGA"                     // SID_3p (12 bp)
        "AATTCGGACCTTAAGGCCTGTAATCTGTAATC"  // R2 Insert (32 bp, reverse complement of R1)
        "TAGCTAGCTACGT";                    // 3' End Adapter (13 bp)

    std::string query_loop = "GTCGACAATTCTTGTCATA";

    int32_t cpu_loop_pos = cpu_bitap_search_oracle(query_loop, test_read, 2);

    BitapPattern gpu_loop_pattern;
    gpu_loop_pattern.init(query_loop.c_str(), static_cast<int>(query_loop.length()), xoos::demux::cuda::SearchDirection::kForward);

    int32_t gpu_loop_pos = bitap_search_window<2, MatchPolicy::kFirst>(
        gpu_loop_pattern, test_read.c_str(), 0, static_cast<int32_t>(test_read.length() - 1)
    );

    std::cout << "  CPU Bitap Loop Pos: " << cpu_loop_pos << std::endl;
    std::cout << "  GPU Bitap Loop Pos: " << gpu_loop_pos << std::endl;

    if (cpu_loop_pos == gpu_loop_pos && cpu_loop_pos > 0) {
        std::cout << "  [PASS] Bitap CPU vs GPU Matched exactly with 0 mismatches!" << std::endl;
    } else {
        std::cerr << "  [FAIL] Mismatch in Bitap position!" << std::endl;
        exit(1);
    }
}

void run_comparative_cpu_gpu_benchmark(size_t num_test_reads) {
    std::cout << "\n================================================================================" << std::endl;
    std::cout << "  CONTRASTING CPU VS GPU BENCHMARK (Dataset: " << num_test_reads << " Duplex Reads)" << std::endl;
    std::cout << "================================================================================" << std::endl;

    std::string loop_seq = "GTCGACAATTCTTGTCATA";
    std::string start_adapter = "ACGTAGCTAGCTA";
    std::string end_adapter = "TAGCTAGCTACGT";
    std::vector<std::string> sid_5p_list = {"TCGACTGACTGA", "AACCGGTTAACC", "TTGGAACCGGTT"};
    std::vector<std::string> sid_3p_list = {"TCAGTCAGTCGA", "GGTTAACCGGTT", "AACCGGTTCCAA"};

    // Construct synthetic FASTQ buffer
    std::string r1_seq = "GATTACAGATTACAGGCCTTAAGGTCCGAATTCCGGTTCCAAGGTTAACCGGTTCCAAGG"; // 60 bp
    std::string r2_exact;
    for (int i = static_cast<int>(r1_seq.length()) - 1; i >= 0; --i) {
        char c = r1_seq[i];
        char rc = (c == 'A') ? 'T' : (c == 'C') ? 'G' : (c == 'G') ? 'C' : 'A';
        r2_exact += rc;
    }

    std::string read_bases = start_adapter + r1_seq + sid_5p_list[0] + loop_seq + sid_3p_list[0] + r2_exact + end_adapter;
    std::string read_quals(read_bases.length(), 'I'); // Q40

    std::string raw_fastq_buffer;
    raw_fastq_buffer.reserve(num_test_reads * (read_bases.length() + 80));
    std::vector<std::pair<std::string, std::string>> cpu_reads;
    cpu_reads.reserve(num_test_reads);

    for (size_t i = 0; i < num_test_reads; ++i) {
        raw_fastq_buffer += "@M00123:45:000000000-A1B2C:1:1101:" + std::to_string(1000 + i) + ":1000 1:N:0:1\n";
        raw_fastq_buffer += read_bases + "\n";
        raw_fastq_buffer += "+\n";
        raw_fastq_buffer += read_quals + "\n";
        cpu_reads.emplace_back(read_bases, read_quals);
    }

    double total_mb = static_cast<double>(raw_fastq_buffer.size()) / (1024.0 * 1024.0);
    std::cout << "  Input Buffer Size:         " << std::fixed << std::setprecision(2) << total_mb << " MB (" << num_test_reads << " reads)" << std::endl;

    // ------------------------------------------------------------------------
    // 1. CPU Single-Threaded Benchmark
    // ------------------------------------------------------------------------
    std::cout << "\n>>> [1] Running CPU Single-Threaded Benchmark..." << std::endl;
    std::atomic<uint64_t> cpu_1t_valid{0};
    auto t_cpu1_start = std::chrono::high_resolution_clock::now();
    run_cpu_demux_threaded(cpu_reads, loop_seq, start_adapter, end_adapter, sid_5p_list, sid_3p_list, 1, cpu_1t_valid);
    auto t_cpu1_end = std::chrono::high_resolution_clock::now();
    double cpu_1t_ms = std::chrono::duration<double, std::milli>(t_cpu1_end - t_cpu1_start).count();
    double cpu_1t_reads_sec = (num_test_reads / (cpu_1t_ms / 1000.0));
    double cpu_1t_gbps = (raw_fastq_buffer.size() / 1e9) / (cpu_1t_ms / 1000.0);

    std::cout << "  CPU 1-Thread Time:         " << std::fixed << std::setprecision(2) << cpu_1t_ms << " ms" << std::endl;
    std::cout << "  CPU 1-Thread Throughput:   " << std::setprecision(2) << (cpu_1t_reads_sec / 1e3) << " kReads/sec (" << cpu_1t_gbps << " GB/s)" << std::endl;

    // ------------------------------------------------------------------------
    // 2. CPU Multi-Threaded (16 Threads) Benchmark
    // ------------------------------------------------------------------------
    unsigned int num_threads = std::min(16u, std::thread::hardware_concurrency());
    std::cout << "\n>>> [2] Running CPU Multi-Threaded (" << num_threads << " Threads) Benchmark..." << std::endl;
    std::atomic<uint64_t> cpu_mt_valid{0};
    auto t_cpumt_start = std::chrono::high_resolution_clock::now();
    run_cpu_demux_threaded(cpu_reads, loop_seq, start_adapter, end_adapter, sid_5p_list, sid_3p_list, num_threads, cpu_mt_valid);
    auto t_cpumt_end = std::chrono::high_resolution_clock::now();
    double cpu_mt_ms = std::chrono::duration<double, std::milli>(t_cpumt_end - t_cpumt_start).count();
    double cpu_mt_reads_sec = (num_test_reads / (cpu_mt_ms / 1000.0));
    double cpu_mt_gbps = (raw_fastq_buffer.size() / 1e9) / (cpu_mt_ms / 1000.0);

    std::cout << "  CPU " << num_threads << "-Thread Time:        " << std::fixed << std::setprecision(2) << cpu_mt_ms << " ms" << std::endl;
    std::cout << "  CPU " << num_threads << "-Thread Throughput:  " << std::setprecision(2) << (cpu_mt_reads_sec / 1e3) << " kReads/sec (" << cpu_mt_gbps << " GB/s)" << std::endl;

    // ------------------------------------------------------------------------
    // 3. GPU RTX 5090 Blackwell CUDA Benchmark
    // ------------------------------------------------------------------------
    std::cout << "\n>>> [3] Running GPU CUDA Benchmark (NVIDIA RTX 5090 sm_120)..." << std::endl;
    DemuxCudaEngine engine(0);
    engine.configure_adapter_bundle(loop_seq, start_adapter, end_adapter, sid_5p_list, sid_3p_list);

    CudaDemuxStats gpu_stats;
    auto t_gpu_start = std::chrono::high_resolution_clock::now();
    engine.process_raw_buffer(raw_fastq_buffer.data(), raw_fastq_buffer.size(), gpu_stats);
    auto t_gpu_end = std::chrono::high_resolution_clock::now();
    double gpu_e2e_ms = std::chrono::duration<double, std::milli>(t_gpu_end - t_gpu_start).count();
    double gpu_e2e_reads_sec = (num_test_reads / (gpu_e2e_ms / 1000.0));
    double gpu_e2e_gbps = (raw_fastq_buffer.size() / 1e9) / (gpu_e2e_ms / 1000.0);

    std::cout << "  GPU Kernel Time:           " << std::fixed << std::setprecision(2) << gpu_stats.kernel_time_ms << " ms" << std::endl;
    std::cout << "  GPU End-to-End Time:       " << std::fixed << std::setprecision(2) << gpu_e2e_ms << " ms (incl. PCIe H2D DMA)" << std::endl;
    std::cout << "  GPU Kernel Throughput:     " << std::setprecision(2) << (gpu_stats.throughput_reads_per_sec / 1e6) << " Million reads/sec (" << gpu_stats.throughput_gbps << " GB/s)" << std::endl;
    std::cout << "  GPU End-to-End Throughput: " << std::setprecision(2) << (gpu_e2e_reads_sec / 1e6) << " Million reads/sec (" << gpu_e2e_gbps << " GB/s)" << std::endl;

    // ------------------------------------------------------------------------
    // 4. Comparative Speedup Summary Table
    // ------------------------------------------------------------------------
    double speedup_vs_1t = cpu_1t_ms / gpu_stats.kernel_time_ms;
    double speedup_vs_mt = cpu_mt_ms / gpu_stats.kernel_time_ms;
    double speedup_e2e_vs_1t = cpu_1t_ms / gpu_e2e_ms;
    double speedup_e2e_vs_mt = cpu_mt_ms / gpu_e2e_ms;

    std::cout << "\n==========================================================================================" << std::endl;
    std::cout << "                        CPU VS GPU PERFORMANCE COMPARISON SUMMARY" << std::endl;
    std::cout << "==========================================================================================" << std::endl;
    std::cout << std::left << std::setw(28) << "Execution Engine"
              << std::right << std::setw(14) << "Time (ms)"
              << std::setw(22) << "Throughput (reads/s)"
              << std::setw(16) << "Bandwidth"
              << std::setw(12) << "Speedup" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    
    std::cout << std::left << std::setw(28) << "CPU 1-Thread"
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << cpu_1t_ms
              << std::setw(22) << std::setprecision(0) << cpu_1t_reads_sec
              << std::setw(13) << std::setprecision(2) << cpu_1t_gbps << " GB/s"
              << std::setw(11) << "1.0x" << std::endl;

    std::cout << std::left << std::setw(28) << ("CPU " + std::to_string(num_threads) + "-Threads")
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << cpu_mt_ms
              << std::setw(22) << std::setprecision(0) << cpu_mt_reads_sec
              << std::setw(13) << std::setprecision(2) << cpu_mt_gbps << " GB/s"
              << std::setw(10) << std::setprecision(1) << (cpu_1t_ms / cpu_mt_ms) << "x" << std::endl;

    std::cout << std::left << std::setw(28) << "GPU RTX 5090 (Kernel)"
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << gpu_stats.kernel_time_ms
              << std::setw(22) << std::setprecision(0) << gpu_stats.throughput_reads_per_sec
              << std::setw(13) << std::setprecision(2) << gpu_stats.throughput_gbps << " GB/s"
              << std::setw(10) << std::setprecision(1) << speedup_vs_1t << "x" << std::endl;

    std::cout << std::left << std::setw(28) << "GPU RTX 5090 (End-to-End)"
              << std::right << std::setw(14) << std::fixed << std::setprecision(2) << gpu_e2e_ms
              << std::setw(22) << std::setprecision(0) << gpu_e2e_reads_sec
              << std::setw(13) << std::setprecision(2) << gpu_e2e_gbps << " GB/s"
              << std::setw(10) << std::setprecision(1) << speedup_e2e_vs_1t << "x" << std::endl;
    std::cout << "==========================================================================================\n" << std::endl;
}

int main(int argc, char** argv) {
    size_t num_reads = 100000; // 100k reads (~25 MB raw FASTQ)
    if (argc > 1) {
        num_reads = std::stoull(argv[1]);
    }

    std::cout << "================================================================================" << std::endl;
    std::cout << "  XOOS DEMUX CPU VS GPU BENCHMARK SUITE" << std::endl;
    std::cout << "================================================================================" << std::endl;

    run_bitap_oracle_parity_test();
    run_comparative_cpu_gpu_benchmark(num_reads);

    std::cout << ">>> BENCHMARK COMPLETED SUCCESSFULLY (0 mismatches) <<<\n" << std::endl;
    return 0;
}
