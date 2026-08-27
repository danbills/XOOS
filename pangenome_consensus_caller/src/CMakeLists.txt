# ============================================================================
# XOOS Pangenome Consensus Caller Core Library CMakeLists.txt
# ============================================================================

find_package(xoos-htslib QUIET)
if (xoos-htslib_FOUND)
    add_default_library(pangenome_consensus_caller_core
        pangenome-consensus-caller-core.cpp
    )
    target_include_directories(pangenome_consensus_caller_core SYSTEM INTERFACE .)
    target_include_directories(pangenome_consensus_caller_core PRIVATE .)
    target_link_libraries(pangenome_consensus_caller_core
        PRIVATE
        xoos-htslib::xoos-htslib
    )
endif ()

# CUDA GPU Acceleration for pangenome_consensus_caller
include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS "CUDA Acceleration: ENABLED for pangenome_consensus_caller_core (Compiler: ${CMAKE_CUDA_COMPILER})")

    add_default_library(pangenome_consensus_caller_cuda_core
        cuda/pangenome_consensus_caller_cuda.cu
    )
    target_include_directories(pangenome_consensus_caller_cuda_core PRIVATE . cuda)
    target_include_directories(pangenome_consensus_caller_cuda_core SYSTEM INTERFACE . cuda)
    target_link_libraries(pangenome_consensus_caller_cuda_core PUBLIC CUDA::cudart)
    target_link_xoos_libraries(pangenome_consensus_caller_cuda_core PUBLIC types)
    target_compile_definitions(pangenome_consensus_caller_cuda_core PUBLIC XOOS_CUDA_ENABLED=1)
endif ()
