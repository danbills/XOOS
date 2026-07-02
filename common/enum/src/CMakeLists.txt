add_default_library(xoos_enum
        xoos/enum/enum-util.cpp
)
target_include_directories(xoos_enum
        PRIVATE
        .
)
target_include_directories(xoos_enum
        SYSTEM
        INTERFACE
        .
)

find_package(fmt REQUIRED)
find_package(magic_enum REQUIRED)

target_link_libraries(xoos_enum
        PUBLIC
        fmt::fmt
        magic_enum::magic_enum
)
target_link_xoos_libraries(xoos_enum
        PUBLIC
        log
)
