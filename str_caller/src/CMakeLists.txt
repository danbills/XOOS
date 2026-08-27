# ============================================================================
# XOOS Short Tandem Repeat (STR) Caller Core Library CMakeLists.txt
# ============================================================================

include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS "CUDA Acceleration: ENABLED for str_caller_core (Compiler: ${CMAKE_CUDA_COMPILER})")

    add_default_library(str_caller_cuda_core
        cuda/str_caller_cuda.cu
    )
    target_include_directories(str_caller_cuda_core PRIVATE . cuda)
    target_include_directories(str_caller_cuda_core SYSTEM INTERFACE . cuda)
    target_link_libraries(str_caller_cuda_core PUBLIC CUDA::cudart)
    target_link_xoos_libraries(str_caller_cuda_core PUBLIC types)
    target_compile_definitions(str_caller_cuda_core PUBLIC XOOS_CUDA_ENABLED=1)
endif ()
