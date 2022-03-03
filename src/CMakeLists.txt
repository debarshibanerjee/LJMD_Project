cmake_minimum_required(VERSION 3.16)

project(ljmd VERSION 1.0 LANGUAGES C)

add_compile_definitions(LJMD_VERSION=1.0)
find_package(MPI REQUIRED)

if(ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
endif()

add_library(ljmd_lib ${CMAKE_SOURCE_DIR}/src/utilities.c ${CMAKE_SOURCE_DIR}/src/verlet.c ${CMAKE_SOURCE_DIR}/src/input.c ${CMAKE_SOURCE_DIR}/src/output.c ${CMAKE_SOURCE_DIR}/src/force.c ${CMAKE_SOURCE_DIR}/src/energy.c)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
link_directories(${CMAKE_SOURCE_DIR}/src)
add_executable(ljmd.x  ${CMAKE_SOURCE_DIR}/src/main.c)

if(ENABLE_OPENMP)
	target_link_libraries(ljmd.x PRIVATE ${MPI_C_LIBRARIES} ljmd_lib m OpenMP::OpenMP_C)
else()
	target_link_libraries(ljmd.x PRIVATE ${MPI_C_LIBRARIES} ljmd_lib m)
endif()

target_compile_options(ljmd.x PRIVATE -O3 -ffast-math -fexpensive-optimizations -msse3 -Wall -Wextra)
