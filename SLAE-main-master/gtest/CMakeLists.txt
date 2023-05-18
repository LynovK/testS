# 'lib' is the folder with Google Test sources
cmake_minimum_required(VERSION 3.0.0)
project(slae)

set(CMAKE_CXX_STANDART 17)

find_package(GTest REQUIRED)

add_executable(Google_Tests_run ThreeDiagTest.cpp CompressedMatrixTest.cpp IterationsTest.cpp MPI_Chebishev.cpp GausZeid_mpi_YacobiTest.cpp "(r,v)_to_Kepler_test.cpp" fastest_grad_descent_test.cpp Related_directions_method.cpp trash.cpp KR2.cpp)

add_test(Google_Tests_run COMMAND ThreeDiagTest.cpp)

#include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
#
## Передаем исполняемый файл с тестами из папки тестов
#add_executable(Google_Tests_run ThreeDiagTest.cpp)

target_link_libraries(Google_Tests_run GTest::gtest_main)