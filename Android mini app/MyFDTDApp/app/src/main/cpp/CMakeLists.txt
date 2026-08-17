cmake_minimum_required(VERSION 3.10)

project(fdtd)

add_library(fdtd SHARED
        ${CMAKE_CURRENT_SOURCE_DIR}/fdtd.cpp
)

find_library(log-lib log)

target_link_libraries(fdtd ${log-lib})

