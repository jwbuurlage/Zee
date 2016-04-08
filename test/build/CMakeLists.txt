cmake_minimum_required(VERSION 3.0.0)

project(zee_test)

include_directories(
    "../../include"
    "../../ext/jwutil/include"
    "../../ext/catch/include"
)

set(TEST_NAME "zee_test")

set(
    TEST_SOURCES
    "../catch.cpp"
    "../linear_algebra.cpp"
    "../loading_matrices.cpp"
)

set(
    LIB_NAMES
    "pthread"
)

add_executable(${TEST_NAME} ${TEST_SOURCES})
target_link_libraries( ${TEST_NAME} ${LIB_NAMES} )
