cmake_minimum_required(VERSION 3.16)

project(ljmd VERSION 1.0 LANGUAGES C)

add_compile_definitions(LJMD_VERSION=1.0)
find_package(MPI)
add_library(testlib ${CMAKE_SOURCE_DIR}/src/utilities.c ${CMAKE_SOURCE_DIR}/src/verlet.c ${CMAKE_SOURCE_DIR}/src/input.c ${CMAKE_SOURCE_DIR}/src/output.c ${CMAKE_SOURCE_DIR}/src/force.c ${CMAKE_SOURCE_DIR}/src/energy.c)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
link_directories(${CMAKE_SOURCE_DIR}/src)
add_executable(ljmd-serial  ${CMAKE_SOURCE_DIR}/src/main.c)
target_link_libraries(ljmd-serial PRIVATE ${MPI_C_LIBRARIES} testlib m)
target_compile_options(ljmd-serial PRIVATE -O3 -ffast-math -fexpensive-optimizations -msse3 -Wall -Wextra)
