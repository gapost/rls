cmake_minimum_required(VERSION 3.5)

project(rlstst LANGUAGES CXX)

find_package(Armadillo REQUIRED)

#set(ARMADILLO_INCLUDE_DIR C:/armadillo-11.2.1/include)
include_directories(${ARMADILLO_INCLUDE_DIR})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(rlstst RLS_testing.cpp)

target_link_libraries(rlstst ${ARMADILLO_LIBRARIES})