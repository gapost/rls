cmake_minimum_required(VERSION 3.5)

project(rlstst LANGUAGES CXX)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(../source)

add_executable(rlsfilter rlsfilter.cpp)

target_link_libraries(rlsfilter ${ARMADILLO_LIBRARIES})

