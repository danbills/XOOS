# ============================================================================
# XOOS Pangenome Consensus Caller Executables CMakeLists.txt
# ============================================================================

if (TARGET pangenome_consensus_caller_core)
    add_default_executable(pangenome_consensus_caller
        pangenome-consensus-caller.cpp
    )
    target_link_libraries(pangenome_consensus_caller
        PRIVATE
        pangenome_consensus_caller_core
    )
endif ()

if (CMAKE_CUDA_COMPILER)
    add_default_executable(pangenome_consensus_caller_cuda_bench
        pangenome-consensus-caller-cuda-bench.cpp
    )
    target_link_libraries(pangenome_consensus_caller_cuda_bench
        PRIVATE
        pangenome_consensus_caller_cuda_core
    )
endif ()
