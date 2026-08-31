#include "alignment_metrics_cuda.cuh"
#include "coverage_pileup_cuda.cuh"
#include "accuracy_metrics_cuda.cuh"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cmath>

namespace xoos::alignment_metrics::cuda {

AlignmentMetricsCudaEngine::AlignmentMetricsCudaEngine(int device_id) : device_id_(device_id) {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaSetDevice(device_id_ < count ? device_id_ : 0);
    }
    cudaStreamCreate(&stream_);
}

AlignmentMetricsCudaEngine::~AlignmentMetricsCudaEngine() {
    free_workspace();
    if (stream_) { cudaStreamDestroy(stream_); stream_ = nullptr; }
}

void AlignmentMetricsCudaEngine::allocate_workspace(size_t max_reads, size_t buffer_size, size_t max_ref_len) {
    free_workspace();
    max_batch_reads_ = max_reads;
    max_buffer_bytes_ = buffer_size;
    ref_length_ = max_ref_len;

    cudaMalloc(&d_reads_, max_batch_reads_ * sizeof(AlignedReadRecord));
    cudaMalloc(&d_seq_buffer_, max_buffer_bytes_);
    cudaMalloc(&d_offsets_, max_batch_reads_ * sizeof(int32_t));
    cudaMalloc(&d_lens_, max_batch_reads_ * sizeof(int32_t));

    cudaMalloc(&d_ref_genome_, ref_length_ * sizeof(char));
    cudaMalloc(&d_diff_array_, (ref_length_ + 1) * sizeof(int32_t));
    cudaMalloc(&d_depth_array_, ref_length_ * sizeof(uint32_t));

    uint32_t scan_blocks = (ref_length_ + kScanBlockSize - 1) / kScanBlockSize;
    cudaMalloc(&d_block_sums_, (scan_blocks + 1) * sizeof(int32_t));

    cudaMalloc(&d_coverage_hist_, kMaxCoverageDepth * sizeof(uint64_t));
    cudaMalloc(&d_hp_total_, kMaxHomopolymerLen * sizeof(uint64_t));
    cudaMalloc(&d_hp_ins_, kMaxHomopolymerLen * sizeof(uint64_t));
    cudaMalloc(&d_hp_del_, kMaxHomopolymerLen * sizeof(uint64_t));
    cudaMalloc(&d_hp_sub_, kMaxHomopolymerLen * sizeof(uint64_t));
}

void AlignmentMetricsCudaEngine::free_workspace() {
    if (d_reads_) { cudaFree(d_reads_); d_reads_ = nullptr; }
    if (d_seq_buffer_) { cudaFree(d_seq_buffer_); d_seq_buffer_ = nullptr; }
    if (d_offsets_) { cudaFree(d_offsets_); d_offsets_ = nullptr; }
    if (d_lens_) { cudaFree(d_lens_); d_lens_ = nullptr; }
    if (d_ref_genome_) { cudaFree(d_ref_genome_); d_ref_genome_ = nullptr; }
    if (d_diff_array_) { cudaFree(d_diff_array_); d_diff_array_ = nullptr; }
    if (d_depth_array_) { cudaFree(d_depth_array_); d_depth_array_ = nullptr; }
    if (d_block_sums_) { cudaFree(d_block_sums_); d_block_sums_ = nullptr; }
    if (d_coverage_hist_) { cudaFree(d_coverage_hist_); d_coverage_hist_ = nullptr; }
    if (d_hp_total_) { cudaFree(d_hp_total_); d_hp_total_ = nullptr; }
    if (d_hp_ins_) { cudaFree(d_hp_ins_); d_hp_ins_ = nullptr; }
    if (d_hp_del_) { cudaFree(d_hp_del_); d_hp_del_ = nullptr; }
    if (d_hp_sub_) { cudaFree(d_hp_sub_); d_hp_sub_ = nullptr; }
}

bool AlignmentMetricsCudaEngine::load_reference_genome(const std::string& ref_sequence) {
    if (ref_sequence.empty()) return false;
    allocate_workspace(max_batch_reads_, max_buffer_bytes_, ref_sequence.length());
    cudaMemcpyAsync(d_ref_genome_, ref_sequence.data(), ref_sequence.length(), cudaMemcpyHostToDevice, stream_);
    cudaStreamSynchronize(stream_);
    return true;
}

bool AlignmentMetricsCudaEngine::compute_metrics(
    const std::vector<AlignedReadRecord>& h_reads,
    const std::vector<std::string>& h_sequences,
    CoverageSummaryMetrics& out_coverage,
    HpAccuracyMetrics& out_hp,
    ReadAlignmentStats& out_stats,
    MetricsExecutionStats& out_exec_stats
) {
    uint64_t num_reads = h_reads.size();
    if (num_reads == 0 || ref_length_ == 0) return false;

    auto t_start = std::chrono::high_resolution_clock::now();

    // 1. Flatten sequence buffer
    std::string flat_seqs;
    std::vector<int32_t> offsets(num_reads);
    std::vector<int32_t> lengths(num_reads);
    size_t cur_off = 0;

    for (uint64_t i = 0; i < num_reads; ++i) {
        offsets[i] = static_cast<int32_t>(cur_off);
        lengths[i] = static_cast<int32_t>(h_sequences[i].length());
        flat_seqs += h_sequences[i];
        cur_off += h_sequences[i].length();
    }

    // 2. Transfer Reads and Clear GPU Accumulators
    cudaMemcpyAsync(d_reads_, h_reads.data(), num_reads * sizeof(AlignedReadRecord), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_seq_buffer_, flat_seqs.data(), flat_seqs.size(), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_offsets_, offsets.data(), num_reads * sizeof(int32_t), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_lens_, lengths.data(), num_reads * sizeof(int32_t), cudaMemcpyHostToDevice, stream_);

    cudaMemsetAsync(d_diff_array_, 0, (ref_length_ + 1) * sizeof(int32_t), stream_);
    cudaMemsetAsync(d_depth_array_, 0, ref_length_ * sizeof(uint32_t), stream_);
    cudaMemsetAsync(d_coverage_hist_, 0, kMaxCoverageDepth * sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_hp_total_, 0, kMaxHomopolymerLen * sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_hp_ins_, 0, kMaxHomopolymerLen * sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_hp_del_, 0, kMaxHomopolymerLen * sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_hp_sub_, 0, kMaxHomopolymerLen * sizeof(uint64_t), stream_);

    // GPU Timing
    cudaEvent_t k_start, k_stop;
    cudaEventCreate(&k_start);
    cudaEventCreate(&k_stop);
    cudaEventRecord(k_start, stream_);

    // 3. Difference Array Deposit Kernel
    int threads = 256;
    int blocks_reads = (num_reads + threads - 1) / threads;
    if (blocks_reads > 170 * 8) blocks_reads = 170 * 8;
    if (blocks_reads == 0) blocks_reads = 1;

    xoos_deposit_difference_array_kernel<<<blocks_reads, threads, 0, stream_>>>(
        d_reads_, num_reads, ref_length_, d_diff_array_
    );

    // 4. Native 2-Pass Parallel Prefix Scan across Reference Genome
    uint32_t scan_blocks = (ref_length_ + kScanBlockSize - 1) / kScanBlockSize;
    xoos_scan_local_blocks_kernel<<<scan_blocks, kScanBlockSize, 0, stream_>>>(
        d_diff_array_, ref_length_, d_depth_array_, d_block_sums_
    );
    xoos_scan_block_sums_kernel<<<1, 1, 0, stream_>>>(
        d_block_sums_, scan_blocks
    );
    xoos_add_block_sums_kernel<<<scan_blocks, kScanBlockSize, 0, stream_>>>(
        d_block_sums_, ref_length_, d_depth_array_
    );

    // 5. Coverage Histogram & Reduction Kernel
    uint64_t* d_total_aligned_bases = nullptr;
    uint64_t* d_covered_bases = nullptr;
    uint32_t* d_max_coverage = nullptr;
    cudaMalloc(&d_total_aligned_bases, sizeof(uint64_t));
    cudaMalloc(&d_covered_bases, sizeof(uint64_t));
    cudaMalloc(&d_max_coverage, sizeof(uint32_t));
    cudaMemsetAsync(d_total_aligned_bases, 0, sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_covered_bases, 0, sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_max_coverage, 0, sizeof(uint32_t), stream_);

    int blocks_ref = (ref_length_ + threads - 1) / threads;
    if (blocks_ref > 170 * 8) blocks_ref = 170 * 8;
    if (blocks_ref == 0) blocks_ref = 1;

    xoos_coverage_histogram_kernel<<<blocks_ref, threads, 0, stream_>>>(
        d_depth_array_, ref_length_, d_coverage_hist_,
        d_total_aligned_bases, d_covered_bases, d_max_coverage
    );

    // 6. Homopolymer Accuracy Kernel
    xoos_hp_accuracy_kernel<<<blocks_reads, threads, 0, stream_>>>(
        d_reads_, d_seq_buffer_, d_offsets_, d_lens_,
        d_ref_genome_, num_reads, ref_length_,
        d_hp_total_, d_hp_ins_, d_hp_del_, d_hp_sub_
    );

    // 7. Global Read Stats Kernel
    uint64_t *d_mapped, *d_unmapped, *d_dup, *d_fwd, *d_rev, *d_mapq_sum, *d_len_sum;
    cudaMalloc(&d_mapped, sizeof(uint64_t));
    cudaMalloc(&d_unmapped, sizeof(uint64_t));
    cudaMalloc(&d_dup, sizeof(uint64_t));
    cudaMalloc(&d_fwd, sizeof(uint64_t));
    cudaMalloc(&d_rev, sizeof(uint64_t));
    cudaMalloc(&d_mapq_sum, sizeof(uint64_t));
    cudaMalloc(&d_len_sum, sizeof(uint64_t));

    cudaMemsetAsync(d_mapped, 0, sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_unmapped, 0, sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_dup, 0, sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_fwd, 0, sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_rev, 0, sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_mapq_sum, 0, sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_len_sum, 0, sizeof(uint64_t), stream_);

    xoos_read_stats_kernel<<<blocks_reads, threads, 0, stream_>>>(
        d_reads_, num_reads, d_mapped, d_unmapped, d_dup, d_fwd, d_rev, d_mapq_sum, d_len_sum
    );

    cudaEventRecord(k_stop, stream_);
    cudaEventSynchronize(k_stop);

    float gpu_kernel_time_ms = 0.0f;
    cudaEventElapsedTime(&gpu_kernel_time_ms, k_start, k_stop);
    cudaEventDestroy(k_start);
    cudaEventDestroy(k_stop);

    // 8. Download Results
    uint64_t h_total_aligned = 0, h_covered = 0;
    uint32_t h_max_cov = 0;
    cudaMemcpy(&h_total_aligned, d_total_aligned_bases, sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_covered, d_covered_bases, sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_max_cov, d_max_coverage, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(out_coverage.coverage_histogram, d_coverage_hist_, kMaxCoverageDepth * sizeof(uint64_t), cudaMemcpyDeviceToHost);

    cudaMemcpy(out_hp.hp_total_bases, d_hp_total_, kMaxHomopolymerLen * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(out_hp.hp_insertion_errors, d_hp_ins_, kMaxHomopolymerLen * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(out_hp.hp_deletion_errors, d_hp_del_, kMaxHomopolymerLen * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(out_hp.hp_substitution_errors, d_hp_sub_, kMaxHomopolymerLen * sizeof(uint64_t), cudaMemcpyDeviceToHost);

    uint64_t h_mapped = 0, h_unmapped = 0, h_dup = 0, h_fwd = 0, h_rev = 0, h_mapq_sum = 0, h_len_sum = 0;
    cudaMemcpy(&h_mapped, d_mapped, sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_unmapped, d_unmapped, sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_dup, d_dup, sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_fwd, d_fwd, sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_rev, d_rev, sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_mapq_sum, d_mapq_sum, sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_len_sum, d_len_sum, sizeof(uint64_t), cudaMemcpyDeviceToHost);

    cudaFree(d_total_aligned_bases);
    cudaFree(d_covered_bases);
    cudaFree(d_max_coverage);
    cudaFree(d_mapped); cudaFree(d_unmapped); cudaFree(d_dup);
    cudaFree(d_fwd); cudaFree(d_rev); cudaFree(d_mapq_sum); cudaFree(d_len_sum);

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    // 9. Compute Statistical Summaries
    out_coverage.total_reference_bases = ref_length_;
    out_coverage.covered_bases = h_covered;
    out_coverage.total_aligned_bases = h_total_aligned;
    out_coverage.mean_coverage = static_cast<double>(h_total_aligned) / ref_length_;
    out_coverage.max_coverage = h_max_cov;
    out_coverage.pct_bases_ge_1x = (static_cast<double>(h_covered) / ref_length_) * 100.0;

    for (size_t i = 0; i < kMaxHomopolymerLen; ++i) {
        if (out_hp.hp_total_bases[i] > 0) {
            uint64_t err_sum = out_hp.hp_insertion_errors[i] + out_hp.hp_deletion_errors[i] + out_hp.hp_substitution_errors[i];
            out_hp.hp_error_rate[i] = (static_cast<double>(err_sum) / out_hp.hp_total_bases[i]) * 100.0;
        }
    }

    out_stats.total_reads = num_reads;
    out_stats.mapped_reads = h_mapped;
    out_stats.unmapped_reads = h_unmapped;
    out_stats.duplicate_reads = h_dup;
    out_stats.forward_strand_reads = h_fwd;
    out_stats.reverse_strand_reads = h_rev;
    out_stats.mean_mapq = (h_mapped > 0) ? (static_cast<double>(h_mapq_sum) / h_mapped) : 0.0;
    out_stats.mean_read_length = (h_mapped > 0) ? (static_cast<double>(h_len_sum) / h_mapped) : 0.0;

    out_exec_stats.total_reads_processed = num_reads;
    out_exec_stats.total_reference_bases_scanned = ref_length_;
    out_exec_stats.total_time_ms = gpu_kernel_time_ms;
    if (gpu_kernel_time_ms > 0) {
        out_exec_stats.throughput_reads_per_sec = num_reads / (gpu_kernel_time_ms / 1000.0);
        out_exec_stats.throughput_gbps = (flat_seqs.size() / 1e9) / (gpu_kernel_time_ms / 1000.0);
    }

    return true;
}

void AlignmentMetricsCudaEngine::print_device_info() const {
    int count = 0;
    cudaGetDeviceCount(&count);
    if (count > 0) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, device_id_ < count ? device_id_ : 0);
        std::cout << "\n================================================================================" << std::endl;
        std::cout << "  XOOS ALIGNMENT METRICS CUDA ENGINE (sm_" << prop.major << prop.minor << ")" << std::endl;
        std::cout << "================================================================================" << std::endl;
        std::cout << "  Device Name:             " << prop.name << std::endl;
        std::cout << "  Compute Capability:      " << prop.major << "." << prop.minor << std::endl;
        std::cout << "  Streaming Multiprocessors: " << prop.multiProcessorCount << " SMs" << std::endl;
        std::cout << "  Total Global Memory:     " << (prop.totalGlobalMem >> 20) << " MB" << std::endl;
        std::cout << "  L2 Cache Size:           " << (prop.l2CacheSize >> 20) << " MB" << std::endl;
        std::cout << "================================================================================\n" << std::endl;
    }
}

} // namespace xoos::alignment_metrics::cuda
