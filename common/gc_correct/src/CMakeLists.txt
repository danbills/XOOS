find_package(Armadillo REQUIRED)
find_package(Taskflow REQUIRED)
add_default_library(xoos_gc_correct
	xoos/gc_correct/gc-correct-taskflow-graph.cpp
	xoos/gc_correct/loess-helper-functions.cpp
	xoos/gc_correct/loess-taskflow-graph.cpp
	xoos/gc_correct/loess.cpp
	xoos/gc_correct/gc-correct.cpp
)
target_include_directories(xoos_gc_correct
        PRIVATE
        .
)
target_include_directories(xoos_gc_correct
        SYSTEM
        INTERFACE
        .
)

target_link_libraries(xoos_gc_correct
        PUBLIC
	armadillo::armadillo
	Taskflow::Taskflow
)

target_link_xoos_libraries(xoos_gc_correct
        PRIVATE
        log
)
