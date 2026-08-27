# ============================================================================
# XOOS Fused Super-Kernel Engine Core Library CMakeLists.txt
# ============================================================================

include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS "CUDA & NVRTC Acceleration: ENABLED for fused_core (Compiler: ${CMAKE_CUDA_COMPILER})")

    add_default_library(fused_cuda_core
        ../../common/fused/fused_stage2_engine.cu
        ../../common/fused/nvrtc_jit_engine.cpp
    )
    target_include_directories(fused_cuda_core PRIVATE ../../common ../../common/fused .)
    target_include_directories(fused_cuda_core SYSTEM INTERFACE ../../common ../../common/fused .)
    target_link_libraries(fused_cuda_core
        PUBLIC
        CUDA::cudart
        CUDA::nvrtc
        ${CMAKE_DL_LIBS}
    )
    target_link_xoos_libraries(fused_cuda_core PUBLIC types)
    target_compile_definitions(fused_cuda_core PUBLIC XOOS_CUDA_ENABLED=1)
endif ()
