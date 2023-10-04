function (printvar var)
    message("${var}: ${${var}}")
endfunction()

cmake_minimum_required(VERSION 3.10)
project(Fast-BVE VERSION 1.0 LANGUAGES CXX C)
set(default_build_type "Release")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (BUILD_TYPE EQUAL Release)
  add_compile_options(-O0)
else ()
  add_compile_options(-g)
endif()

printvar(BUILD_TYPE)

find_package(MPI REQUIRED COMPONENTS C)
find_package(LAPACK REQUIRED)

add_subdirectory(src)
add_subdirectory(executables)
