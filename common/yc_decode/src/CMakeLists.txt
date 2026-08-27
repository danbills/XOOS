find_package(xoos-htslib QUIET)

if (xoos-htslib_FOUND)
    add_default_library(xoos_yc_decode
        xoos/yc-decode/yc-decoder.cpp
        xoos/yc-decode/cigar-utils.cpp
    )

    target_include_directories(xoos_yc_decode PRIVATE .)
    target_include_directories(xoos_yc_decode SYSTEM INTERFACE .)

    target_link_libraries(xoos_yc_decode
        PUBLIC
        xoos-htslib::xoos-htslib
        PRIVATE
    )

    target_link_xoos_libraries(xoos_yc_decode
        PUBLIC
        error
        io
        log
        types
        PRIVATE
        enum
        util
    )
endif ()
