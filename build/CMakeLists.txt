cmake_minimum_required(VERSION 3.0.0)

project(zee)

include_directories(
    "../include"
    "../ext/jwutil/include"
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wfatal-errors -Wextra -g -O3")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

set(
    EXAMPLE_DIR
    "../examples/"
)

set(
    EXAMPLE_SOURCES
    "gmres_experiments.cpp"
    "linear_algebra.cpp"
    "partitioner.cpp"
    "pulp_experiments.cpp"
    "vector_partitioner.cpp"
)

set(
    LIB_NAMES
    "pthread"
)

foreach(source_file ${EXAMPLE_SOURCES})
    string(REPLACE ".cpp" "" source_name ${source_file})
    add_executable(${source_name} ${EXAMPLE_DIR}${source_file})
    target_link_libraries( ${source_name} ${LIB_NAMES} )
endforeach(source_file)

add_subdirectory(../test/build ../test/build)
