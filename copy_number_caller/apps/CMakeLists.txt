# ============================================================================
# XOOS Copy Number Caller Executables CMakeLists.txt
# ============================================================================

if (TARGET copy_number_caller_core)
    add_default_executable(copy_number_caller
        copy-number-caller.cpp
    )
    target_link_libraries(copy_number_caller
        PRIVATE
        copy_number_caller_core
    )
    target_link_xoos_libraries(copy_number_caller
        PRIVATE
        enum
        log
        util
        cli
        io
    )

    add_default_executable(predict_sex
        predict-sex.cpp
    )
    target_link_libraries(predict_sex
        PRIVATE
        copy_number_caller_core
    )
    target_link_xoos_libraries(predict_sex
        PRIVATE
        enum
        log
        util
        cli
        io
        types
    )
endif ()

if (TARGET generate_panel_of_normals_core)
    add_default_executable(generate_panel_of_normals
        generate-panel-of-normals.cpp
    )
    target_link_libraries(generate_panel_of_normals
        PRIVATE
        generate_panel_of_normals_core
    )
    target_link_xoos_libraries(generate_panel_of_normals
        PRIVATE
        enum
        log
        util
        cli
        io
        types
    )
endif ()

if (CMAKE_CUDA_COMPILER)
    add_default_executable(copy_number_caller_cuda_bench
        copy-number-caller-cuda-bench.cpp
    )
    target_link_libraries(copy_number_caller_cuda_bench
        PRIVATE
        copy_number_caller_cuda_core
    )
endif ()
