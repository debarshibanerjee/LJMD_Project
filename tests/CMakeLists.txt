
cmake_minimum_required(VERSION 3.16)
project(GTESTDEMO VERSION 1.0 LANGUAGES C)
include_directories(${CMAKE_SOURCE_DIR}/include)
find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
add_definitions(-DOMPI_SKIP_MPICXX)
option(ENABLE_TESTING "Enable building unit tests" OFF)
if(ENABLE_TESTING)
  enable_testing()
  enable_language(CXX)

  include(${CMAKE_SOURCE_DIR}/tests/GTest.cmake)
  add_executable(test_verlet test_verlet.cpp)
  target_include_directories(test_verlet PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_verlet PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME Verlet COMMAND ${MPIEXEC_EXECUTABLE} -np 1 test_verlet)

  add_executable(test_utils test_utils.cpp)
  target_include_directories(test_utils PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(test_utils PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME Utils COMMAND test_utils)

  add_executable(test_energy test_energy.cpp)
  target_include_directories(test_energy PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_energy PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME Energy COMMAND ${MPIEXEC_EXECUTABLE} -np 1 test_energy)

  add_executable(test_ekin test_ekin.cpp)
  target_include_directories(test_ekin PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_ekin PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME Ekin COMMAND ${MPIEXEC_EXECUTABLE} -np 1 test_ekin)

  add_executable(test_verlet_integration_detailed test_verlet_integration_detailed.cpp)
  target_include_directories(test_verlet_integration_detailed PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_verlet_integration_detailed PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME Detailed_Verlet_Integration_Test COMMAND ${MPIEXEC_EXECUTABLE} -np 1 test_verlet_integration_detailed)

  add_executable(test_pbc test_pbc.cpp)
  target_include_directories(test_pbc PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_pbc PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME PBC COMMAND test_pbc)

  add_executable(test_input test_input.cpp)
  target_include_directories(test_input PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_input PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME Input COMMAND test_input)

  add_executable(test_output test_output.cpp)
  target_include_directories(test_output PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_output PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME output COMMAND test_output)

  add_executable(test_mpi test_mpi.cpp)
  target_include_directories(test_mpi PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_mpi PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME mpi COMMAND ${MPIEXEC_EXECUTABLE} --oversubscribe -np 4 test_mpi)

  add_executable(test_force test_force.cpp)
  target_include_directories(test_force PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME Force COMMAND ${MPIEXEC_EXECUTABLE} --oversubscribe -np 4 test_force)

  add_executable(test_force_3LawN test_force_3LawN.cpp)
  target_include_directories(test_force_3LawN PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_force_3LawN PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME Force_3LawN_1st_version COMMAND ${MPIEXEC_EXECUTABLE} -np 1 test_force_3LawN)

  add_executable(test_force_3LawN_more_opt test_force_3LawN_more_opt.cpp)
  target_include_directories(test_force_3LawN_more_opt PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_force_3LawN_more_opt PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME Force_3LawN_more_opt_2nd_version COMMAND ${MPIEXEC_EXECUTABLE} -np 1 test_force_3LawN_more_opt)

  add_executable(test_allocate_arrays test_allocate_arrays.cpp)
  target_include_directories(test_allocate_arrays PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_allocate_arrays PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME Array_Allocation COMMAND test_allocate_arrays)

  add_executable(test_free_arrays test_free_arrays.cpp)
  target_include_directories(test_free_arrays PRIVATE ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
  target_link_libraries(test_free_arrays PRIVATE GTest::GTestMain GTest::GTest ${MPI_C_LIBRARIES} testlib)
  add_test(NAME Array_Freeing COMMAND test_free_arrays)
endif()

