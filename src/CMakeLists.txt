cmake_minimum_required(VERSION 3.16)

project(ljmd VERSION 1.0 LANGUAGES C)

add_compile_definitions(LJMD_VERSION=1.0)

if(ENABLE_SIMPLE_FORCE)
	add_definitions(-DSIMPLE)
endif()

find_package(MPI REQUIRED)

if(ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
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
