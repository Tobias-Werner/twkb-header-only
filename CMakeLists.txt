PROJECT(TWKB-Header-Only)

cmake_minimum_required(VERSION 3.9.0)
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra -Wconversion")


add_executable(TWKB-Tests test/src/TWKB_Tests.cpp test/include/catch.hpp include/TWKB.h)

target_include_directories(TWKB-Tests PRIVATE include test/include)



