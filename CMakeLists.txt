function (printvar var)
    message("${var}: ${${var}}")
endfunction()

cmake_minimum_required(VERSION 3.10)
project(FastBVE VERSION 1.0 LANGUAGES CXX C)
set(DEFAULT_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
option(BUILD_SHARED_LIBS "Build shared libraries, rather than static" ON)

if (BUILD_SHARED_LIBS)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()
set(CMAKE_CXX_FLAGS_RELEASE -O3)

# turn off deprecated env variable handling:
# cmake will now use $<Package>_ROOT env variables in its search for Package
# dependencies
CMAKE_POLICY(SET CMP0074 NEW)

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.24)
  # turn off deprecated cmake download timestamp handling
  # cmake will now set the local timestamp of a dependency to match its package
  # creation timestamp, rather than its download timestamp.
  CMAKE_POLICY(SET CMP0135 NEW)
endif()
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

# CMake auto-generated config settings will go in this file
configure_file(fastbve-config.h.in fastbve-config.h)

# The following commands prepare FastBVE to be packaged as a library and used
# by other codes as a dependency.
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# generate config file for other CMake projects to find this project
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/FastBVEConfig.cmake
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FastBVE"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
# include version info so that client codes can ensure compatibility
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/FastBVEConfigVersion.cmake"
  VERSION "${FastBVE_VERSION_MAJOR}.${FastBVE_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion
)
# install configured package files
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/FastBVEConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/FastBVEConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FastBVE
)
# export configured package files to the build directory, too; this allows
# clients to use FastBVE from a FastBVE build directory (without having to
# install it).
#export(EXPORT FastBVETargets
#  FILE "${CMAKE_CURRENT_BINARY_DIR}/FastBVETargets.cmake"
#)
#install(EXPORT FastBVETargets
#  FILE FastBVETargets.cmake
#  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FastBVE
#)
