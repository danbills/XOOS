# ============================================================================
# XOOS CUDA Build Configuration Module
# File: common/cmake/setup_cuda.cmake
# ============================================================================

macro(setup_cuda)
    option(ENABLE_CUDA "Enable NVIDIA CUDA GPU acceleration" ON)
    if (ENABLE_CUDA)
        include(CheckLanguage)
        check_language(CUDA)
        if (CMAKE_CUDA_COMPILER)
            enable_language(CUDA)
            set(CMAKE_CUDA_STANDARD 17)
            set(CMAKE_CUDA_STANDARD_REQUIRED ON)
            find_package(CUDAToolkit REQUIRED)
            set(XOOS_CUDA_AVAILABLE TRUE)
            message(STATUS "CUDA Acceleration: ENABLED (Compiler: ${CMAKE_CUDA_COMPILER})")
            add_compile_definitions(XOOS_CUDA_ENABLED=1)
        else ()
            set(XOOS_CUDA_AVAILABLE FALSE)
            message(STATUS "CUDA Acceleration: DISABLED (No CUDA compiler found; CPU-only)")
        endif ()
    else ()
        set(XOOS_CUDA_AVAILABLE FALSE)
        message(STATUS "CUDA Acceleration: DISABLED by user (ENABLE_CUDA=OFF)")
    endif ()
endmacro()
