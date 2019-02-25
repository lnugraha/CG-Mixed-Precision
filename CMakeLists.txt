cmake_minimum_required(VERSION 3.12)
project(HACS)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(${DCMTK_INCLUDE_DIRS})
include_directories(include)

set(SOURCE_FILES main.cpp)
file(GLOB LAPLACE "src/laplace.cpp")
file(GLOB LIBMATH "src/libmath.cpp")

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

add_executable(HACS ${SOURCE_FILES} ${LAPLACE} ${LIBMATH})
