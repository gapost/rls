cmake_minimum_required(VERSION 3.5)

project(test LANGUAGES CXX)

find_package(Armadillo REQUIRED)

include_directories(${ARMADILLO_INCLUDE_DIR})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(../source)

add_executable(filtertest filtertest.cpp)

target_link_libraries(filtertest ${ARMADILLO_LIBRARIES})
