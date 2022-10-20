cmake_minimum_required(VERSION 3.5)

project(test LANGUAGES CXX)

find_package(Eigen3 3.3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(../source)

add_executable(chol_testing chol_testing.cpp)

target_link_libraries(chol_testing Eigen3::Eigen)
