function (printvar var)
    message("${var}: ${${var}}")
endfunction()

cmake_minimum_required(VERSION 3.10)
project(Fast-BVE VERSION 1.0 LANGUAGES CXX C)
set(default_build_type "Release")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

printvar(PROJECT_SOURCE_DIR)

find_package(MPI REQUIRED COMPONENTS C)
find_package(LAPACK REQUIRED)

add_subdirectory(src)
add_subdirectory(executables)
