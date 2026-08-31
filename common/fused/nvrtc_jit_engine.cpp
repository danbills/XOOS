#include "nvrtc_jit_engine.hpp"
#include <nvrtc.h>
#include <cuda.h>
#include <dlfcn.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstdlib>
#include <sys/stat.h>
#include <iomanip>

namespace xoos::fused::jit {

typedef CUresult (*pfn_cuInit)(unsigned int);
typedef CUresult (*pfn_cuModuleLoadData)(CUmodule*, const void*);
typedef CUresult (*pfn_cuModuleGetFunction)(CUfunction*, CUmodule, const char*);
typedef CUresult (*pfn_cuLaunchKernel)(CUfunction, unsigned int, unsigned int, unsigned int,
                                       unsigned int, unsigned int, unsigned int,
                                       unsigned int, CUstream, void**, void**);

static void* s_cuda_handle = nullptr;
static pfn_cuInit s_cuInit = nullptr;
static pfn_cuModuleLoadData s_cuModuleLoadData = nullptr;
static pfn_cuModuleGetFunction s_cuModuleGetFunction = nullptr;
static pfn_cuLaunchKernel s_cuLaunchKernel = nullptr;

static bool load_cuda_driver() {
    if (s_cuda_handle) return true;

    const char* lib_names[] = {
        "libcuda.so.1",
        "libcuda.so",
        "/usr/lib/wsl/lib/libcuda.so.1",
        "/usr/lib64/libcuda.so.1",
        "/usr/lib/x86_64-linux-gnu/libcuda.so.1",
        "/lib/libcuda.so.1",
        "/lib/stubs/libcuda.so"
    };

    for (const char* name : lib_names) {
        s_cuda_handle = dlopen(name, RTLD_LAZY | RTLD_GLOBAL);
        if (s_cuda_handle) break;
    }

    if (!s_cuda_handle) return false;

    s_cuInit = reinterpret_cast<pfn_cuInit>(dlsym(s_cuda_handle, "cuInit"));
    s_cuModuleLoadData = reinterpret_cast<pfn_cuModuleLoadData>(dlsym(s_cuda_handle, "cuModuleLoadData"));
    s_cuModuleGetFunction = reinterpret_cast<pfn_cuModuleGetFunction>(dlsym(s_cuda_handle, "cuModuleGetFunction"));
    s_cuLaunchKernel = reinterpret_cast<pfn_cuLaunchKernel>(dlsym(s_cuda_handle, "cuLaunchKernel"));

    return (s_cuInit && s_cuModuleLoadData && s_cuModuleGetFunction && s_cuLaunchKernel);
}

NvrtcSuperKernelJit::NvrtcSuperKernelJit() {
    if (load_cuda_driver() && s_cuInit) {
        CUresult res = s_cuInit(0);
        driver_initialized_ = (res == CUDA_SUCCESS);
    }

    const char* home = std::getenv("HOME");
    if (home) {
        cache_dir_ = std::string(home) + "/.cache/xoos/kernels";
    } else {
        cache_dir_ = "/tmp/xoos_kernels";
    }
    ensure_cache_dir();
}

NvrtcSuperKernelJit::~NvrtcSuperKernelJit() {}

void NvrtcSuperKernelJit::ensure_cache_dir() {
    mkdir((std::string(std::getenv("HOME") ? std::getenv("HOME") : "/tmp") + "/.cache").c_str(), 0777);
    mkdir((std::string(std::getenv("HOME") ? std::getenv("HOME") : "/tmp") + "/.cache/xoos").c_str(), 0777);
    mkdir(cache_dir_.c_str(), 0777);
}

std::string NvrtcSuperKernelJit::compute_hash(const DynamicPolicyConfig& config) const {
    std::ostringstream ss;
    ss << "k_fused_s2_smem_"
       << (config.enable_rescue ? "1" : "0") << "_"
       << (config.enable_collapse ? "1" : "0") << "_"
       << (config.enable_gc_metrics ? "1" : "0") << "_"
       << (config.enable_insert_metrics ? "1" : "0") << "_"
       << config.min_family_size << "_"
       << static_cast<int>(config.adjusted_bq);
    return ss.str();
}

std::string NvrtcSuperKernelJit::generate_source(const DynamicPolicyConfig& config) const {
    std::ostringstream src;
    src << R"(
typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int uint32_t;
typedef unsigned long long uint64_t;

#define MAX_SEQ_LEN 256
#define NUM_GC_BINS 101
#define MAX_INSERT_SIZE 1000

struct __align__(16) FusedReadRecord {
    uint64_t read_id;
    uint64_t barcode_hash;
    uint32_t chr_id;
    uint32_t pos;
    uint16_t length;
    uint16_t insert_size;
    uint8_t mapq;
    uint8_t is_reverse_strand;
    uint8_t has_yc_tag;
    uint8_t is_duplicate;
    float r1_graph_score;
    float r2_graph_score;
    char sequence[MAX_SEQ_LEN];
    char yc_tag[MAX_SEQ_LEN];
    uint8_t base_qual[MAX_SEQ_LEN];
};

struct __align__(16) GlobalMetricsAccumulator {
    uint64_t total_reads;
    uint64_t total_bases;
    uint64_t total_gc_bases;
    uint64_t total_rescued_reads;
    uint64_t total_base_corrections;
    uint64_t total_collapsed_families;
    uint64_t total_duplicates_marked;
    uint32_t gc_histogram[NUM_GC_BINS];
    uint32_t insert_size_histogram[MAX_INSERT_SIZE];
};

__device__ __forceinline__ char decode_yc(char yc, char r1_base) {
    switch (yc) {
        case 'c': return (r1_base == 'A') ? 'C' : (r1_base == 'C' ? 'A' : r1_base);
        case 'g': return (r1_base == 'A') ? 'G' : (r1_base == 'G' ? 'A' : r1_base);
        case 't': return (r1_base == 'A') ? 'T' : (r1_base == 'T' ? 'A' : r1_base);
        case 'k': return (r1_base == 'C') ? 'G' : (r1_base == 'G' ? 'C' : r1_base);
        case 'y': return (r1_base == 'C') ? 'T' : (r1_base == 'T' ? 'C' : r1_base);
        case 'w': return (r1_base == 'G') ? 'T' : (r1_base == 'T' ? 'G' : r1_base);
        case 'C': return 'C';
        case 'G': return 'G';
        case 'T': return 'T';
        case 'A': return 'A';
        default: return r1_base;
    }
}

extern "C" __global__ void xoos_dynamic_fused_stage2_kernel(
    FusedReadRecord* __restrict__ d_reads,
    uint64_t num_reads,
    GlobalMetricsAccumulator* __restrict__ d_metrics
) {
    __shared__ uint32_t s_gc_hist[NUM_GC_BINS];
    __shared__ unsigned long long s_total_bases;
    __shared__ unsigned long long s_total_gc_bases;
    __shared__ unsigned long long s_rescued_reads;
    __shared__ unsigned long long s_base_corrections;
    __shared__ unsigned long long s_collapsed_families;
    __shared__ unsigned long long s_duplicates_marked;

    uint32_t tid = threadIdx.x;
    if (tid < NUM_GC_BINS) s_gc_hist[tid] = 0;
    if (tid == 0) {
        s_total_bases = 0;
        s_total_gc_bases = 0;
        s_rescued_reads = 0;
        s_base_corrections = 0;
        s_collapsed_families = 0;
        s_duplicates_marked = 0;
    }
    __syncthreads();

    size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t stride = (size_t)gridDim.x * blockDim.x;

    for (size_t r = idx; r < num_reads; r += stride) {
        FusedReadRecord& read = d_reads[r];
        uint16_t len = read.length;
        if (len > MAX_SEQ_LEN) len = MAX_SEQ_LEN;

)";

    if (config.enable_rescue) {
        src << "        // 1. Dynamic JIT Pangenome Rescue\n"
            << "        if (read.has_yc_tag && (read.r2_graph_score > read.r1_graph_score)) {\n"
            << "            uint32_t corr = 0;\n"
            << "            for (uint16_t i = 0; i < len; ++i) {\n"
            << "                char yc = read.yc_tag[i];\n"
            << "                if (yc != '*' && yc != '~' && yc != '\\0') {\n"
            << "                    char r2 = decode_yc(yc, read.sequence[i]);\n"
            << "                    if (r2 != read.sequence[i]) {\n"
            << "                        read.sequence[i] = r2;\n"
            << "                        read.base_qual[i] = " << static_cast<int>(config.adjusted_bq) << ";\n"
            << "                        corr++;\n"
            << "                    }\n"
            << "                }\n"
            << "            }\n"
            << "            if (corr > 0) {\n"
            << "                atomicAdd(&s_rescued_reads, 1ULL);\n"
            << "                atomicAdd(&s_base_corrections, (unsigned long long)corr);\n"
            << "            }\n"
            << "        }\n";
    }

    if (config.enable_gc_metrics) {
        src << "        // 2. Dynamic JIT GC Metrics\n"
            << "        uint32_t gc_cnt = 0;\n"
            << "        for (uint16_t i = 0; i < len; ++i) {\n"
            << "            char b = read.sequence[i];\n"
            << "            if (b == 'G' || b == 'C' || b == 'g' || b == 'c') gc_cnt++;\n"
            << "        }\n"
            << "        uint32_t gc_pct = (len > 0) ? (gc_cnt * 100 / len) : 0;\n"
            << "        if (gc_pct > 100) gc_pct = 100;\n"
            << "        atomicAdd(&s_gc_hist[gc_pct], 1);\n"
            << "        atomicAdd(&s_total_gc_bases, (unsigned long long)gc_cnt);\n"
            << "        atomicAdd(&s_total_bases, (unsigned long long)len);\n";
    }

    if (config.enable_insert_metrics) {
        src << "        // 3. Dynamic JIT Insert Size Metrics\n"
            << "        if (read.insert_size < MAX_INSERT_SIZE) {\n"
            << "            atomicAdd(&d_metrics->insert_size_histogram[read.insert_size], 1);\n"
            << "        }\n";
    }

    if (config.enable_collapse) {
        src << "        // 4. Dynamic JIT Barcode Collapsing\n"
            << "        if (read.barcode_hash != 0) {\n"
            << "            if (read.is_duplicate) {\n"
            << "                atomicAdd(&s_duplicates_marked, 1ULL);\n"
            << "            } else {\n"
            << "                atomicAdd(&s_collapsed_families, 1ULL);\n"
            << "            }\n"
            << "        }\n";
    }

    src << R"(
    }

    __syncthreads();

    if (tid < NUM_GC_BINS && s_gc_hist[tid] > 0) {
        atomicAdd(&d_metrics->gc_histogram[tid], s_gc_hist[tid]);
    }

    if (tid == 0) {
        if (s_total_bases > 0) {
            atomicAdd((unsigned long long*)&d_metrics->total_bases, s_total_bases);
            atomicAdd((unsigned long long*)&d_metrics->total_gc_bases, s_total_gc_bases);
        }
        if (s_rescued_reads > 0) {
            atomicAdd((unsigned long long*)&d_metrics->total_rescued_reads, s_rescued_reads);
            atomicAdd((unsigned long long*)&d_metrics->total_base_corrections, s_base_corrections);
        }
        if (s_collapsed_families > 0) {
            atomicAdd((unsigned long long*)&d_metrics->total_collapsed_families, s_collapsed_families);
        }
        if (s_duplicates_marked > 0) {
            atomicAdd((unsigned long long*)&d_metrics->total_duplicates_marked, s_duplicates_marked);
        }
    }
}
)";
    return src.str();
}

bool NvrtcSuperKernelJit::get_or_compile_kernel(
    const DynamicPolicyConfig& config,
    CUmodule& out_module,
    CUfunction& out_function,
    double& out_compile_time_ms
) {
    if (!load_cuda_driver() || !s_cuModuleLoadData || !s_cuModuleGetFunction) {
        return false;
    }

    std::string kernel_hash = compute_hash(config);
    std::string ptx_path = cache_dir_ + "/" + kernel_hash + ".ptx";

    // 1. Check disk cache
    std::ifstream cached_file(ptx_path, std::ios::binary);
    if (cached_file.is_open()) {
        std::stringstream ss;
        ss << cached_file.rdbuf();
        std::string ptx_data = ss.str();

        CUresult res = s_cuModuleLoadData(&out_module, ptx_data.c_str());
        if (res == CUDA_SUCCESS) {
            res = s_cuModuleGetFunction(&out_function, out_module, "xoos_dynamic_fused_stage2_kernel");
            if (res == CUDA_SUCCESS) {
                out_compile_time_ms = 0.0; // Cached hit!
                return true;
            }
        }
    }

    // 2. Synthesize & Compile in memory via NVRTC
    auto t_start = std::chrono::high_resolution_clock::now();

    std::string src = generate_source(config);
    nvrtcProgram prog;
    nvrtcResult nres = nvrtcCreateProgram(&prog, src.c_str(), "xoos_fused.cu", 0, nullptr, nullptr);
    if (nres != NVRTC_SUCCESS) return false;

    const char* opts[] = {
        "--std=c++17",
        "--gpu-architecture=compute_80",
        "--use_fast_math"
    };

    nres = nvrtcCompileProgram(prog, 3, opts);
    if (nres != NVRTC_SUCCESS) {
        size_t log_size;
        nvrtcGetProgramLogSize(prog, &log_size);
        std::vector<char> log(log_size);
        nvrtcGetProgramLog(prog, log.data());
        std::cerr << ">>> NVRTC JIT Compilation Error:\n" << log.data() << std::endl;
        nvrtcDestroyProgram(&prog);
        return false;
    }

    size_t ptx_size;
    nvrtcGetPTXSize(prog, &ptx_size);
    std::vector<char> ptx(ptx_size);
    nvrtcGetPTX(prog, ptx.data());
    nvrtcDestroyProgram(&prog);

    auto t_end = std::chrono::high_resolution_clock::now();
    out_compile_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    // 3. Save to disk cache
    std::ofstream out_cache(ptx_path, std::ios::binary);
    if (out_cache.is_open()) {
        out_cache.write(ptx.data(), ptx.size());
    }

    // 4. Load module into driver
    CUresult res = s_cuModuleLoadData(&out_module, ptx.data());
    if (res != CUDA_SUCCESS) return false;

    res = s_cuModuleGetFunction(&out_function, out_module, "xoos_dynamic_fused_stage2_kernel");
    return (res == CUDA_SUCCESS);
}

bool NvrtcSuperKernelJit::launch_kernel(
    CUfunction func,
    unsigned int gridDimX,
    unsigned int blockDimX,
    CUstream stream,
    void** args
) {
    if (!load_cuda_driver() || !s_cuLaunchKernel) return false;

    CUresult res = s_cuLaunchKernel(
        func,
        gridDimX, 1, 1,
        blockDimX, 1, 1,
        0,
        stream,
        args,
        nullptr
    );
    return (res == CUDA_SUCCESS);
}

} // namespace xoos::fused::jit
