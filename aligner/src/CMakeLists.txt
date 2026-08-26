# ============================================================================
# XOOS Aligner Core Library CMakeLists.txt
# ============================================================================

include(CheckLanguage)
check_language(CUDA)

if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    message(STATUS "CUDA Acceleration: ENABLED for aligner_core (Compiler: ${CMAKE_CUDA_COMPILER})")

    add_default_library(aligner_core
        cuda/aligner_cuda.cu
    )

    target_include_directories(aligner_core PRIVATE . cuda)
    target_include_directories(aligner_core SYSTEM INTERFACE . cuda)

    target_link_xoos_libraries(aligner_core
        PUBLIC
        types
        util
        error
        log
    )
    target_compile_definitions(aligner_core PUBLIC XOOS_CUDA_ENABLED=1)
endif ()
