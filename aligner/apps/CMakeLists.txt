# ============================================================================
# XOOS Aligner Executables CMakeLists.txt
# ============================================================================

if (CMAKE_CUDA_COMPILER)
    add_default_executable(aligner_cuda_bench
        aligner-cuda-bench.cpp
    )
    target_link_libraries(aligner_cuda_bench
        PRIVATE
        aligner_core
    )
endif ()

