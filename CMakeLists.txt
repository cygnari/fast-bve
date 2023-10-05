function (printvar var)
    message("${var}: ${${var}}")
endfunction()

cmake_minimum_required(VERSION 3.10)
project(Fast-BVE VERSION 1.0 LANGUAGES CXX C)
set(default_build_type Release)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS -O3)

printvar(BUILD_TYPE)

CMAKE_POLICY(SET CMP0135 NEW)
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
FetchContent_MakeAvailable(googletest)

find_package(MPI REQUIRED COMPONENTS C)
find_package(LAPACK REQUIRED)

enable_testing()
add_subdirectory(src)
add_subdirectory(executables)
add_subdirectory(tests)
