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

    find_package(PkgConfig REQUIRED)
    pkg_check_modules(HTSLIB REQUIRED htslib)

    add_default_executable(roche_sbx_cuda_eval
        ../../apps/roche-sbx-cuda-eval.cpp
    )
    target_include_directories(roche_sbx_cuda_eval PRIVATE ../../common ../../common/fused ${HTSLIB_INCLUDE_DIRS})
    target_link_libraries(roche_sbx_cuda_eval
        PRIVATE
        fused_cuda_core
        ${HTSLIB_LIBRARIES}
    )
endif ()
