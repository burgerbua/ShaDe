cmake_minimum_required (VERSION 3.6.2)
project (ShaDe)

include_directories(${CMAKE_SOURCE_DIR}/sources/)

file (GLOB headers "${CMAKE_SOURCE_DIR}/sources/*.hxx")
source_group ("Headers" FILES ${headers})

add_executable (shade_test main.cpp ${headers})