add_default_library(xoos_error
    xoos/error/error.cpp
)
target_include_directories(xoos_error
    PRIVATE
    .
)
target_include_directories(xoos_error
    SYSTEM
    INTERFACE
    .
)

find_package(fmt REQUIRED)

target_link_libraries(xoos_error
    PUBLIC
    fmt::fmt
)
