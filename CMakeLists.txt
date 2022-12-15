cmake_minimum_required(VERSION 3.5)

# set path to additional CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/filtertester ${CMAKE_MODULE_PATH})

project(rls LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory( filter )

# add_subdirectory( filtertester )
