cmake_minimum_required(VERSION 3.10)

project(ljmd VERSION 1.0 LANGUAGES C)

set(exe_name "ljmd-serial.x")


include_directories(${CMAKE_SOURCE_DIR}/include)
link_directories(${CMAKE_SOURCE_DIR}/src)
add_executable(${exe_name} ${CMAKE_SOURCE_DIR}/src/main.c)
add_compile_options(-Wall -g -O3 -ffast-math -fomit-frame-pointer)


