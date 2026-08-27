#include "demux_cuda.cuh"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <stdexcept>

namespace xoos::demux::cuda {

namespace {

/**
 * @brief CUDA Kernel: Warp-Parallel Duplex Demultiplexing and Consensus Collapsing
 *
 * Each thread or warp processes one sequencing read:
 *   1. Scans FASTQ boundaries (@name, seq, +, qual)
 *   2. Evaluates Bitap 64-bit state machine in ALU registers to locate Loop & SIDs
 *   3. Trims Start/End adapters and extracts R1 and R2 spans
 *   4. Performs pairwise R1 <-> R2_RC consensus and Q-score recalculation
 */
__global__ void xoos_demux_duplex_kernel(
    const char* __restrict__ d_buffer,
    size_t buffer_size,
    const GpuAdapterBundle* __restrict__ d_bundle,
    DuplexTrimResult* __restrict__ d_out_trims,
    ConsensusReadResult* __restrict__ d_out_consensus,
    uint32_t* __restrict__ d_global_read_count,
    uint32_t* __restrict__ d_valid_duplex_count,
    uint64_t* __restrict__ d_total_concordant_bases,
    uint64_t* __restrict__ d_total_consensus_bases,
    uint32_t max_reads
) {
    size_t tid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    size_t stride = static_cast<size_t>(gridDim.x) * blockDim.x;

    for (size_t pos = tid; pos < buffer_size; pos += stride) {
        // Fast line synchronization on FASTQ '@' record start
        if (d_buffer[pos] == '@' && (pos == 0 || d_buffer[pos - 1] == '\n')) {
            size_t p = pos;

            // Line 1: Header
            while (p < buffer_size && d_buffer[p] != '\n') p++;
            if (p >= buffer_size) break;
            p++;

            // Line 2: Sequence
            size_t seq_start = p;
            while (p < buffer_size && d_buffer[p] != '\n' && d_buffer[p] != '\r') p++;
            size_t seq_len = p - seq_start;
            if (p >= buffer_size) break;
            while (p < buffer_size && (d_buffer[p] == '\n' || d_buffer[p] == '\r')) p++;

            // Line 3: '+' Separator
            if (p >= buffer_size || d_buffer[p] != '+') continue;
            while (p < buffer_size && d_buffer[p] != '\n') p++;
            if (p >= buffer_size) break;
            p++;

            // Line 4: Qualities
            size_t qual_start = p;
            while (p < buffer_size && d_buffer[p] != '\n' && d_buffer[p] != '\r') p++;
            size_t qual_len = p - qual_start;

            if (seq_len == 0 || seq_len != qual_len || seq_len < static_cast<size_t>(kMinDuplexReadLen)) {
                continue;
            }

            // Warp-level atomic aggregation to eliminate SM atomic memory contention
            uint32_t active_lanes = __activemask();
            int lane_id = threadIdx.x & 31;
            uint32_t total_in_warp = __popc(active_lanes);
            uint32_t warp_base_idx = 0;

            if (lane_id == (__ffs(active_lanes) - 1)) {
                warp_base_idx = atomicAdd(d_global_read_count, total_in_warp);
            }
            warp_base_idx = __shfl_sync(active_lanes, warp_base_idx, __ffs(active_lanes) - 1);
            uint32_t r_idx = warp_base_idx + __popc(active_lanes & ((1U << lane_id) - 1));

            if (r_idx >= max_reads) {
                break;
            }

            const char* seq = d_buffer + seq_start;
            const char* qual = d_buffer + qual_start;

            // 1. In-Register Hairpin & SID Extraction
            DuplexTrimResult trim;
            trim.read_idx = r_idx;
            extract_duplex_hairpin(*d_bundle, seq, static_cast<int32_t>(seq_len), trim);

            // 2. In-Register Duplex Consensus & Quality Recalibration
            ConsensusReadResult cons;
            collapse_duplex_consensus(seq, qual, trim, cons);

            // 3. Write results to output SoA buffers
            d_out_trims[r_idx] = trim;
            d_out_consensus[r_idx] = cons;

            if (trim.is_valid_duplex && cons.consensus_len > 0) {
                atomicAdd(d_valid_duplex_count, 1U);
                atomicAdd(reinterpret_cast<unsigned long long*>(d_total_concordant_bases),
                          static_cast<unsigned long long>(cons.concordant_bases));
                atomicAdd(reinterpret_cast<unsigned long long*>(d_total_consensus_bases),
                          static_cast<unsigned long long>(cons.consensus_len));
            }
        }
    }
}

} // anonymous namespace

DemuxCudaEngine::DemuxCudaEngine(int device_id) : device_id_(device_id) {
    cudaSetDevice(device_id_);
    cudaStreamCreate(&stream_);
    cudaMalloc(&d_bundle_, sizeof(GpuAdapterBundle));
    allocate_workspace(256 * 1024 * 1024, max_batch_reads_); // 256 MB default buffer
}

DemuxCudaEngine::~DemuxCudaEngine() {
    free_workspace();
    if (d_bundle_) {
        cudaFree(d_bundle_);
        d_bundle_ = nullptr;
    }
    if (stream_) {
        cudaStreamDestroy(stream_);
        stream_ = nullptr;
    }
}

void DemuxCudaEngine::allocate_workspace(size_t buffer_capacity, size_t max_reads) {
    d_raw_buffer_capacity_ = buffer_capacity;
    max_batch_reads_ = max_reads;

    cudaMalloc(&d_raw_buffer_, d_raw_buffer_capacity_);
    cudaMalloc(&d_trim_results_, max_batch_reads_ * sizeof(DuplexTrimResult));
    cudaMalloc(&d_consensus_results_, max_batch_reads_ * sizeof(ConsensusReadResult));
}

void DemuxCudaEngine::free_workspace() {
    if (d_raw_buffer_) { cudaFree(d_raw_buffer_); d_raw_buffer_ = nullptr; }
    if (d_trim_results_) { cudaFree(d_trim_results_); d_trim_results_ = nullptr; }
    if (d_consensus_results_) { cudaFree(d_consensus_results_); d_consensus_results_ = nullptr; }
}

bool DemuxCudaEngine::configure_adapter_bundle(
    const std::string& loop_seq,
    const std::string& start_adapter,
    const std::string& end_adapter,
    const std::vector<std::string>& sid_5p_list,
    const std::vector<std::string>& sid_3p_list
) {
    h_bundle_.loop_fw.init(loop_seq.c_str(), static_cast<int>(loop_seq.length()), SearchDirection::kForward);
    h_bundle_.loop_bw.init(loop_seq.c_str(), static_cast<int>(loop_seq.length()), SearchDirection::kBackward);

    h_bundle_.start_adapter_fw.init(start_adapter.c_str(), static_cast<int>(start_adapter.length()), SearchDirection::kForward);
    h_bundle_.end_adapter_fw.init(end_adapter.c_str(), static_cast<int>(end_adapter.length()), SearchDirection::kForward);

    int num_samples = static_cast<int>(sid_5p_list.size());
    if (num_samples > kMaxSidsPerBundle) {
        num_samples = kMaxSidsPerBundle;
    }
    h_bundle_.num_samples = num_samples;

    for (int i = 0; i < num_samples; ++i) {
        h_bundle_.sid_5p_matchers[i].init(
            sid_5p_list[i].c_str(), static_cast<int>(sid_5p_list[i].length()), SearchDirection::kForward
        );
        h_bundle_.sid_3p_matchers[i].init(
            sid_3p_list[i].c_str(), static_cast<int>(sid_3p_list[i].length()), SearchDirection::kForward
        );
    }

    // Upload configured bundle to device
    cudaMemcpyAsync(d_bundle_, &h_bundle_, sizeof(GpuAdapterBundle), cudaMemcpyHostToDevice, stream_);
    cudaStreamSynchronize(stream_);
    return true;
}

bool DemuxCudaEngine::process_raw_buffer(
    const char* h_raw_buffer,
    size_t buffer_size,
    CudaDemuxStats& out_stats
) {
    if (!h_raw_buffer || buffer_size == 0) {
        return false;
    }

    if (buffer_size > d_raw_buffer_capacity_) {
        free_workspace();
        allocate_workspace(buffer_size * 2, max_batch_reads_);
    }

    // Allocate counters on device
    uint32_t* d_global_read_count = nullptr;
    uint32_t* d_valid_duplex_count = nullptr;
    uint64_t* d_total_concordant_bases = nullptr;
    uint64_t* d_total_consensus_bases = nullptr;

    cudaMalloc(&d_global_read_count, sizeof(uint32_t));
    cudaMalloc(&d_valid_duplex_count, sizeof(uint32_t));
    cudaMalloc(&d_total_concordant_bases, sizeof(uint64_t));
    cudaMalloc(&d_total_consensus_bases, sizeof(uint64_t));

    cudaMemsetAsync(d_global_read_count, 0, sizeof(uint32_t), stream_);
    cudaMemsetAsync(d_valid_duplex_count, 0, sizeof(uint32_t), stream_);
    cudaMemsetAsync(d_total_concordant_bases, 0, sizeof(uint64_t), stream_);
    cudaMemsetAsync(d_total_consensus_bases, 0, sizeof(uint64_t), stream_);

    // 1. DMA Transfer raw buffer to VRAM
    auto t_start = std::chrono::high_resolution_clock::now();
    cudaMemcpyAsync(d_raw_buffer_, h_raw_buffer, buffer_size, cudaMemcpyHostToDevice, stream_);

    // 2. Launch CUDA Demux Kernel
    int threads_per_block = 256;
    int num_blocks = 170 * 4; // High SM occupancy for Blackwell sm_120
    
    xoos_demux_duplex_kernel<<<num_blocks, threads_per_block, 0, stream_>>>(
        d_raw_buffer_, buffer_size, d_bundle_,
        d_trim_results_, d_consensus_results_,
        d_global_read_count, d_valid_duplex_count,
        d_total_concordant_bases, d_total_consensus_bases,
        static_cast<uint32_t>(max_batch_reads_)
    );

    cudaStreamSynchronize(stream_);
    auto t_end = std::chrono::high_resolution_clock::now();

    // 3. Read back statistics
    uint32_t h_global_read_count = 0;
    uint32_t h_valid_duplex_count = 0;
    uint64_t h_total_concordant_bases = 0;
    uint64_t h_total_consensus_bases = 0;

    cudaMemcpy(&h_global_read_count, d_global_read_count, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_valid_duplex_count, d_valid_duplex_count, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_total_concordant_bases, d_total_concordant_bases, sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_total_consensus_bases, d_total_consensus_bases, sizeof(uint64_t), cudaMemcpyDeviceToHost);

    cudaFree(d_global_read_count);
    cudaFree(d_valid_duplex_count);
    cudaFree(d_total_concordant_bases);
    cudaFree(d_total_consensus_bases);

    double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    out_stats.kernel_time_ms = elapsed_ms;
    out_stats.wallclock_time_sec = elapsed_ms / 1000.0;
    out_stats.total_reads_processed = h_global_read_count;
    out_stats.valid_duplex_reads = h_valid_duplex_count;
    out_stats.sample_assigned_reads = h_valid_duplex_count;
    out_stats.failed_hairpin_reads = (h_global_read_count > h_valid_duplex_count) ? (h_global_read_count - h_valid_duplex_count) : 0;
    
    if (out_stats.wallclock_time_sec > 0) {
        out_stats.throughput_reads_per_sec = out_stats.total_reads_processed / out_stats.wallclock_time_sec;
        out_stats.throughput_gbps = (buffer_size / 1e9) / out_stats.wallclock_time_sec;
    }
    if (h_total_consensus_bases > 0) {
        out_stats.mean_concordance_rate = static_cast<float>(h_total_concordant_bases) / static_cast<float>(h_total_consensus_bases);
    }

    return true;
}

void DemuxCudaEngine::print_device_info() const {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device_id_);
    std::cout << "\n================================================================================" << std::endl;
    std::cout << "  XOOS DEMUX CUDA ACCELERATION ENGINE (sm_" << prop.major << prop.minor << ")" << std::endl;
    std::cout << "================================================================================" << std::endl;
    std::cout << "  Device Name:             " << prop.name << std::endl;
    std::cout << "  Compute Capability:      " << prop.major << "." << prop.minor << std::endl;
    std::cout << "  Streaming Multiprocessors: " << prop.multiProcessorCount << " SMs" << std::endl;
    std::cout << "  Total Global Memory:     " << (prop.totalGlobalMem >> 20) << " MB" << std::endl;
    std::cout << "  L2 Cache Size:           " << (prop.l2CacheSize >> 20) << " MB" << std::endl;
    std::cout << "================================================================================\n" << std::endl;
}

} // namespace xoos::demux::cuda
