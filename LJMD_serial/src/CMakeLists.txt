cmake_minimum_required(VERSION 3.15)

project(ljmd VERSION 1.0 LANGUAGES C)

add_compile_definitions(LJMD_VERSION=1.0)


include_directories(${CMAKE_SOURCE_DIR}/include)
link_directories(${CMAKE_SOURCE_DIR}/src)

add_library(ljmd_lib ${CMAKE_SOURCE_DIR}/src/utilities.c ${CMAKE_SOURCE_DIR}/src/verlet.c ${CMAKE_SOURCE_DIR}/src/input.c ${CMAKE_SOURCE_DIR}/src/output.c ${CMAKE_SOURCE_DIR}/src/force.c ${CMAKE_SOURCE_DIR}/src/energy.c)

add_executable(ljmd.x  ${CMAKE_SOURCE_DIR}/src/main.c)

target_compile_options(ljmd.x PRIVATE -O3 -ffast-math -fomit-frame-pointer -fexpensive-optimizations -msse3 -Wall )


target_link_libraries(ljmd.x PRIVATE ljmd_lib m)
