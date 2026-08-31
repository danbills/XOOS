# ============================================================================
# XOOS Tumor Fraction Estimator (TFE) Core Library CMakeLists.txt
# ============================================================================

include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS "CUDA Acceleration: ENABLED for tfe_core (Compiler: ${CMAKE_CUDA_COMPILER})")

    add_default_library(tfe_cuda_core
        cuda/tumor_fraction_cuda.cu
    )
    target_include_directories(tfe_cuda_core PRIVATE . cuda)
    target_include_directories(tfe_cuda_core SYSTEM INTERFACE . cuda)
    target_link_libraries(tfe_cuda_core PUBLIC CUDA::cudart)
    target_link_xoos_libraries(tfe_cuda_core PUBLIC types)
    target_compile_definitions(tfe_cuda_core PUBLIC XOOS_CUDA_ENABLED=1)
endif ()
