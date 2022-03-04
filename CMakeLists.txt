# CMake configuration demo

cmake_minimum_required(VERSION 3.16)
project(ljmd)

if(ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(src)
add_subdirectory(tests)
