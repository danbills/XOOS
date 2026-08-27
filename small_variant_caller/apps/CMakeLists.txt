# ============================================================================
# XOOS Small Variant Caller Executables CMakeLists.txt
# ============================================================================

if (TARGET svc_core)
    add_default_executable(compute_bam_features
        compute-bam-features.cpp
    )
    target_link_libraries(compute_bam_features
        PRIVATE
        svc_core
    )
    target_link_xoos_libraries(compute_bam_features
        PRIVATE
        cli
        log
        util
        io
        enum
    )

    add_default_executable(filter_variants
        filter-variants.cpp
    )
    target_link_libraries(filter_variants
        PRIVATE
        svc_core
    )
    target_link_xoos_libraries(filter_variants
        PRIVATE
        cli
        io
        util
        log
        enum
    )

    add_default_executable(compute_vcf_features
        compute-vcf-features.cpp
    )
    target_link_libraries(compute_vcf_features
        PRIVATE
        svc_core
    )
    target_link_xoos_libraries(compute_vcf_features
        PRIVATE
        cli
        log
        util
        io
        enum
    )

    add_default_executable(train_model
        train-model.cpp
    )
    target_link_libraries(train_model
        PRIVATE
        svc_core
    )
    target_link_xoos_libraries(train_model
        PRIVATE
        cli
        util
        io
        log
        enum
    )

    add_default_executable(vcf_to_bed
        vcf-to-bed.cpp
    )
    target_link_libraries(vcf_to_bed
        PRIVATE
        svc_core
    )
    target_link_xoos_libraries(vcf_to_bed
        PRIVATE
        cli
        io
        util
        log
        enum
    )

    add_default_executable(generate_pon
        generate-pon.cpp
    )
    target_link_libraries(generate_pon
        PRIVATE
        svc_core
    )
    target_link_xoos_libraries(generate_pon
        PRIVATE
        cli
        log
        util
        io
        enum
    )
endif ()

if (CMAKE_CUDA_COMPILER)
    add_default_executable(small_variant_caller_cuda_bench
        small-variant-caller-cuda-bench.cpp
    )
    target_link_libraries(small_variant_caller_cuda_bench
        PRIVATE
        svc_cuda_core
    )
endif ()
