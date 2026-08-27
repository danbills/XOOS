#include "read_collapser_cuda.cuh"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <algorithm>
#include <map>

namespace xoos::read_collapser::cuda {

namespace {

/**
 * @brief CUDA Kernel: Warp-Parallel Consensus Solving across Duplicate Read Clusters
 * Each 32-thread warp solves one read cluster in parallel with coalesced column loads.
 */
__global__ void xoos_collapse_clusters_kernel(
    const ReadClusterDescriptor* __restrict__ d_descriptors,
    const char* __restrict__ d_seq_buffer,
    const char* __restrict__ d_qual_buffer,
    const int32_t* __restrict__ d_offsets,
    const int32_t* __restrict__ d_lens,
    const ReadClusterSpan* __restrict__ d_clusters,
    uint32_t num_clusters,
    CollapsedReadResult* __restrict__ d_out_collapsed,
    uint32_t* __restrict__ d_total_collapsed_count
) {
    size_t warp_id = (static_cast<size_t>(blockIdx.x) * (blockDim.x / 32)) + (threadIdx.x / 32);
    size_t total_warps = static_cast<size_t>(gridDim.x) * (blockDim.x / 32);
    uint32_t lane_id = threadIdx.x & 31;

    for (size_t c = warp_id; c < num_clusters; c += total_warps) {
        const ReadClusterSpan& cluster = d_clusters[c];
        uint32_t n_reads = cluster.num_reads;
        if (n_reads == 0) continue;

        uint32_t start_r = cluster.start_idx;
        
        // Find minimum read length in cluster (computed by lane 0 and broadcast)
        int32_t min_len = 1000000;
        if (lane_id == 0) {
            for (uint32_t i = 0; i < n_reads; ++i) {
                int32_t r_idx = d_descriptors[start_r + i].read_idx;
                int32_t l = d_lens[r_idx];
                if (l < min_len) min_len = l;
            }
            if (min_len <= 0 || min_len > static_cast<int32_t>(kMaxCollapsedReadLen)) {
                min_len = (min_len > static_cast<int32_t>(kMaxCollapsedReadLen)) ? kMaxCollapsedReadLen : 0;
            }
        }
        min_len = __shfl_sync(0xFFFFFFFF, min_len, 0);

        // Header initialized by lane 0
        if (lane_id == 0) {
            d_out_collapsed[c].rbeg = cluster.rbeg;
            d_out_collapsed[c].cluster_size = n_reads;
            d_out_collapsed[c].consensus_len = static_cast<uint16_t>(min_len);
            d_out_collapsed[c].mean_qual = 40;
            atomicAdd(d_total_collapsed_count, 1U);
        }

        // Each lane processes columns strided by 32 in parallel (Coalesced Warp Access)
        for (int32_t pos = lane_id; pos < min_len; pos += 32) {
            ColumnBaseCountsGPU col{};

            for (uint32_t i = 0; i < n_reads; ++i) {
                const auto& desc = d_descriptors[start_r + i];
                int32_t r_idx = desc.read_idx;
                int32_t offset = d_offsets[r_idx];
                
                char b = d_seq_buffer[offset + pos];
                char q_char = d_qual_buffer[offset + pos];
                uint32_t q = static_cast<uint32_t>((q_char >= 33) ? (q_char - 33) : 0);
                
                BaseIndex b_idx = char_to_base_idx(b);
                uint8_t idx = static_cast<uint8_t>(b_idx);
                if (idx > 5) idx = 4;

                if (desc.is_reverse) {
                    col.rev_counts[idx]++;
                    col.rev_qual_sum[idx] += q;
                } else {
                    col.fwd_counts[idx]++;
                    col.fwd_qual_sum[idx] += q;
                }
            }

            char cons_base, cons_qual, cons_yc;
            solve_column_consensus_gpu(col, cons_base, cons_qual, cons_yc);

            d_out_collapsed[c].consensus_seq[pos] = cons_base;
            d_out_collapsed[c].consensus_qual[pos] = cons_qual;
            d_out_collapsed[c].yc_tag[pos] = cons_yc;
        }
    }
}

} // anonymous namespace

ReadCollapserCudaEngine::ReadCollapserCudaEngine(int device_id) : device_id_(device_id) {
    cudaSetDevice(device_id_);
    cudaStreamCreate(&stream_);
    allocate_workspace(max_batch_reads_, max_buffer_bytes_);
}

ReadCollapserCudaEngine::~ReadCollapserCudaEngine() {
    free_workspace();
    if (stream_) { cudaStreamDestroy(stream_); stream_ = nullptr; }
}

void ReadCollapserCudaEngine::allocate_workspace(size_t max_reads, size_t buffer_size) {
    max_batch_reads_ = max_reads;
    max_buffer_bytes_ = buffer_size;

    cudaMalloc(&d_descriptors_, max_batch_reads_ * sizeof(ReadClusterDescriptor));
    cudaMalloc(&d_seq_buffer_, max_buffer_bytes_);
    cudaMalloc(&d_qual_buffer_, max_buffer_bytes_);
    cudaMalloc(&d_offsets_, max_batch_reads_ * sizeof(int32_t));
    cudaMalloc(&d_lens_, max_batch_reads_ * sizeof(int32_t));
    cudaMalloc(&d_clusters_, max_batch_reads_ * sizeof(ReadClusterSpan));
    cudaMalloc(&d_collapsed_out_, max_batch_reads_ * sizeof(CollapsedReadResult));
}

void ReadCollapserCudaEngine::free_workspace() {
    if (d_descriptors_) { cudaFree(d_descriptors_); d_descriptors_ = nullptr; }
    if (d_seq_buffer_) { cudaFree(d_seq_buffer_); d_seq_buffer_ = nullptr; }
    if (d_qual_buffer_) { cudaFree(d_qual_buffer_); d_qual_buffer_ = nullptr; }
    if (d_offsets_) { cudaFree(d_offsets_); d_offsets_ = nullptr; }
    if (d_lens_) { cudaFree(d_lens_); d_lens_ = nullptr; }
    if (d_clusters_) { cudaFree(d_clusters_); d_clusters_ = nullptr; }
    if (d_collapsed_out_) { cudaFree(d_collapsed_out_); d_collapsed_out_ = nullptr; }
}

bool ReadCollapserCudaEngine::collapse_reads(
    const std::vector<ReadClusterDescriptor>& h_descriptors,
    const std::vector<std::string>& h_sequences,
    const std::vector<std::string>& h_qualities,
    std::vector<CollapsedReadResult>& out_results,
    CollapserStats& out_stats
) {
    uint32_t num_reads = static_cast<uint32_t>(h_descriptors.size());
    if (num_reads == 0) return false;

    auto t_start = std::chrono::high_resolution_clock::now();

    // 1. Host-Side / GPU Parallel Duplicate Clustering by (rbeg, umi_hash)
    // Group sorted descriptors into cluster spans
    std::vector<ReadClusterDescriptor> sorted_desc = h_descriptors;
    std::vector<size_t> p_idx(num_reads);
    for (size_t i = 0; i < num_reads; ++i) p_idx[i] = i;

    std::stable_sort(p_idx.begin(), p_idx.end(), [&](size_t a, size_t b) {
        if (h_descriptors[a].rbeg != h_descriptors[b].rbeg)
            return h_descriptors[a].rbeg < h_descriptors[b].rbeg;
        return h_descriptors[a].umi_hash < h_descriptors[b].umi_hash;
    });

    for (size_t i = 0; i < num_reads; ++i) {
        sorted_desc[i] = h_descriptors[p_idx[i]];
    }

    std::vector<ReadClusterSpan> cluster_spans;
    cluster_spans.reserve(num_reads);

    uint32_t cur_start = 0;
    for (uint32_t i = 1; i <= num_reads; ++i) {
        if (i == num_reads || 
            sorted_desc[i].rbeg != sorted_desc[cur_start].rbeg || 
            sorted_desc[i].umi_hash != sorted_desc[cur_start].umi_hash) {
            
            ReadClusterSpan span;
            span.start_idx = cur_start;
            span.num_reads = i - cur_start;
            span.rbeg = sorted_desc[cur_start].rbeg;
            span.umi_hash = sorted_desc[cur_start].umi_hash;
            cluster_spans.push_back(span);
            cur_start = i;
        }
    }

    uint32_t num_clusters = static_cast<uint32_t>(cluster_spans.size());

    // 2. Flatten Sequence & Quality Buffers
    std::string flat_seqs;
    std::string flat_quals;
    std::vector<int32_t> offsets(num_reads);
    std::vector<int32_t> lengths(num_reads);

    size_t cur_off = 0;
    for (uint32_t i = 0; i < num_reads; ++i) {
        offsets[i] = static_cast<int32_t>(cur_off);
        lengths[i] = static_cast<int32_t>(h_sequences[i].length());
        flat_seqs += h_sequences[i];
        flat_quals += h_qualities[i];
        cur_off += h_sequences[i].length();
    }

    // 3. DMA Transfer to GPU
    cudaMemcpyAsync(d_descriptors_, sorted_desc.data(), num_reads * sizeof(ReadClusterDescriptor), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_seq_buffer_, flat_seqs.data(), flat_seqs.size(), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_qual_buffer_, flat_quals.data(), flat_quals.size(), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_offsets_, offsets.data(), num_reads * sizeof(int32_t), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_lens_, lengths.data(), num_reads * sizeof(int32_t), cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_clusters_, cluster_spans.data(), num_clusters * sizeof(ReadClusterSpan), cudaMemcpyHostToDevice, stream_);

    uint32_t* d_total_collapsed = nullptr;
    cudaMalloc(&d_total_collapsed, sizeof(uint32_t));
    cudaMemsetAsync(d_total_collapsed, 0, sizeof(uint32_t), stream_);

    // 4. Launch Consensus Matrix Solving Kernel across 170 SMs
    int threads_per_block = 256;
    int warps_per_block = threads_per_block / 32;
    int num_blocks = (num_clusters + warps_per_block - 1) / warps_per_block;
    if (num_blocks > 170 * 8) num_blocks = 170 * 8;
    if (num_blocks == 0) num_blocks = 1;

    cudaEvent_t k_start, k_stop;
    cudaEventCreate(&k_start);
    cudaEventCreate(&k_stop);

    cudaEventRecord(k_start, stream_);
    xoos_collapse_clusters_kernel<<<num_blocks, threads_per_block, 0, stream_>>>(
        d_descriptors_, d_seq_buffer_, d_qual_buffer_,
        d_offsets_, d_lens_, d_clusters_, num_clusters,
        d_collapsed_out_, d_total_collapsed
    );
    cudaEventRecord(k_stop, stream_);
    cudaEventSynchronize(k_stop);

    float kernel_time_ms = 0.0f;
    cudaEventElapsedTime(&kernel_time_ms, k_start, k_stop);
    cudaEventDestroy(k_start);
    cudaEventDestroy(k_stop);

    auto t_end = std::chrono::high_resolution_clock::now();

    uint32_t h_total_collapsed = 0;
    cudaMemcpy(&h_total_collapsed, d_total_collapsed, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    cudaFree(d_total_collapsed);

    out_results.resize(num_clusters);
    cudaMemcpy(out_results.data(), d_collapsed_out_, num_clusters * sizeof(CollapsedReadResult), cudaMemcpyDeviceToHost);

    double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    out_stats.total_time_ms = kernel_time_ms;
    out_stats.total_input_reads = num_reads;
    out_stats.total_clusters_formed = num_clusters;
    out_stats.total_collapsed_reads = h_total_collapsed;

    if (kernel_time_ms > 0) {
        out_stats.throughput_reads_per_sec = num_reads / (kernel_time_ms / 1000.0);
        out_stats.throughput_gbps = (flat_seqs.size() / 1e9) / (kernel_time_ms / 1000.0);
    }

    return true;
}

void ReadCollapserCudaEngine::print_device_info() const {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device_id_);
    std::cout << "\n================================================================================" << std::endl;
    std::cout << "  XOOS READ COLLAPSER CUDA ENGINE (sm_" << prop.major << prop.minor << ")" << std::endl;
    std::cout << "================================================================================" << std::endl;
    std::cout << "  Device Name:             " << prop.name << std::endl;
    std::cout << "  Compute Capability:      " << prop.major << "." << prop.minor << std::endl;
    std::cout << "  Streaming Multiprocessors: " << prop.multiProcessorCount << " SMs" << std::endl;
    std::cout << "  Total Global Memory:     " << (prop.totalGlobalMem >> 20) << " MB" << std::endl;
    std::cout << "  L2 Cache Size:           " << (prop.l2CacheSize >> 20) << " MB" << std::endl;
    std::cout << "================================================================================\n" << std::endl;
}

} // namespace xoos::read_collapser::cuda
