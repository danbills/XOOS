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
    target_link_xoos_libraries(aligner_cuda_bench
        PRIVATE
        util
        enum
        cli
        concurrent
        log
        io
    )
endif ()
