# ============================================================================
# XOOS Short Tandem Repeat (STR) Caller Executables CMakeLists.txt
# ============================================================================

if (CMAKE_CUDA_COMPILER)
    add_default_executable(str_caller_cuda_bench
        str-caller-cuda-bench.cpp
    )
    target_link_libraries(str_caller_cuda_bench
        PRIVATE
        str_caller_cuda_core
    )
endif ()
