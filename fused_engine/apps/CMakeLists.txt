# ============================================================================
# XOOS Fused Super-Kernel Executables CMakeLists.txt
# ============================================================================

if (CMAKE_CUDA_COMPILER)
    add_default_executable(fused_stage2_cuda_bench
        ../../apps/fused-stage2-cuda-bench.cpp
    )
    target_include_directories(fused_stage2_cuda_bench PRIVATE ../../common ../../common/fused)
    target_link_libraries(fused_stage2_cuda_bench
        PRIVATE
        fused_cuda_core
    )
endif ()
