cmake_minimum_required (VERSION 3.6.2)
project (ShaDe)

set (CMAKE_CXX_STANDARD 11)

include_directories(${CMAKE_SOURCE_DIR}/sources/)

file (GLOB headers "${CMAKE_SOURCE_DIR}/sources/*.hxx")
source_group ("Headers" FILES ${headers})

add_executable (shade_test main.cpp ${headers})
