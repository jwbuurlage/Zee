cmake_minimum_required(VERSION 3.0.0)

project(zee_test)

include_directories(
    "../../include"
    "../../ext/jwutil/include"
    "../../ext/catch/include"
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

set(TEST_NAME "zee_test")

set(
    TEST_SOURCES
    "../catch.cpp"
    "../hypergraphs.cpp"
    "../iterative_solvers.cpp"
    "../linear_algebra.cpp"
    "../loading_matrices.cpp"
    "../partitioners.cpp"
    "../sparse_storage.cpp"
)

set(
    LIB_NAMES
    "pthread"
)

add_executable(${TEST_NAME} ${TEST_SOURCES})
target_link_libraries( ${TEST_NAME} ${LIB_NAMES} )
