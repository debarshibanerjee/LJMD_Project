# CMake configuration demo

cmake_minimum_required(VERSION 3.16)
project(ljmd)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(src)
add_subdirectory(tests)
