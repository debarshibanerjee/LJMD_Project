
cmake_minimum_required(VERSION 3.16)
project(GTESTDEMO VERSION 1.0 LANGUAGES C)
include_directories(${CMAKE_SOURCE_DIR}/include)

option(ENABLE_TESTING "Enable building unit tests" OFF)
if(ENABLE_TESTING)
  enable_testing()
  enable_language(CXX)
  
  include(${CMAKE_SOURCE_DIR}/tests/GTest.cmake)
  add_executable(test_verlet test_verlet.cpp)
  target_include_directories(test_verlet PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(test_verlet PRIVATE GTest::GTestMain GTest::GTest testlib)
  add_test(NAME Verlet COMMAND test_verlet)

  add_executable(test_utils test_utils.cpp)  
  target_include_directories(test_utils PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(test_utils PRIVATE GTest::GTestMain GTest::GTest testlib)
  add_test(NAME Utils COMMAND test_utils)

  add_executable(test_energy test_energy.cpp)  
  target_include_directories(test_energy PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(test_energy PRIVATE GTest::GTestMain GTest::GTest testlib)
  add_test(NAME Energy COMMAND test_energy)

endif()
