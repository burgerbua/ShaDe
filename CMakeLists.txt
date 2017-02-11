cmake_minimum_required (VERSION 3.6.2)
project (ShaDe)

set (CMAKE_CXX_STANDARD 11)

set(EIGEN3_INCLUDE_DIR "/opt/local/include/eigen3" CACHE FILEPATH "Set EIGEN3_INCLUDE_DIR to allow cmake to find Eigen3")
set(VTK_DIR "/opt/local/include/vtk-7.1" CACHE FILEPATH "Set VTK_DIR to allow cmake to find VTK")

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR}/headers/ ${EIGEN3_INCLUDE_DIR})

file (GLOB headers "${CMAKE_SOURCE_DIR}/headers/*.hxx")
source_group ("Headers" FILES ${headers})

add_executable (shade_test ${CMAKE_SOURCE_DIR}/tests/main.cpp ${headers})
