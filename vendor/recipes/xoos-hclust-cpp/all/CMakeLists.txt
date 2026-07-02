cmake_minimum_required(VERSION 3.25)

project(hclust-cpp)

add_library(hclust-cpp fastcluster.cpp)

install(TARGETS hclust-cpp DESTINATION lib)
install(FILES fastcluster.h DESTINATION include)
