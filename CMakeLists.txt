cmake_minimum_required(VERSION 3.16)

project(ljmd_serial VERSION 1.0 LANGUAGES C)

set(exe_name "ljmd-serial.x")

add_executable(${exe_name} src/ljmd.c)

add_compile_options(-Wall -g -O3 -ffast-math -fomit-frame-pointer)

add_compile_definitions(LJMD_VERSION=1.0)

target_link_libraries(${exe_name} m)
