cmake_minimum_required(VERSION 3.25)
include(../common/cmake/setup.cmake)

project(pangenome_consensus_caller)
setup_default_app_project(pangenome_consensus_caller LINT_DISABLE)

if (TARGET pangenome_consensus_caller_core AND TARGET xoos_yc_decode)
    target_link_xoos_libraries(pangenome_consensus_caller_core
        PRIVATE
        yc_decode
        io
    )
endif ()
