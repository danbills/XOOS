# ============================================================================
# XOOS Tumor Fraction Estimator (TFE) Executables CMakeLists.txt
# ============================================================================

if (CMAKE_CUDA_COMPILER)
    add_default_executable(tumor_fraction_estimator_cuda_bench
        tumor-fraction-estimator-cuda-bench.cpp
    )
    target_link_libraries(tumor_fraction_estimator_cuda_bench
        PRIVATE
        tfe_cuda_core
    )
endif ()
