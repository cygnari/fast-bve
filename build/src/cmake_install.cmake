# Install script for directory: /Users/cygnari/documents/Research/FastSummation/fast-bve/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/cygnari/documents/Research/FastSummation/fast-bve/include/amr.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/conservation_fixer.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/direct_sum_utils.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/fast_sum_utils.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/general_utils.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/green_funcs.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/init_utils.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/input_utils.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/interp_utils.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/io_utils.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/mpi_utils.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/rhs_utils.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/structs.hpp;/Users/cygnari/documents/Research/FastSummation/fast-bve/include/vorticity_functions.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/cygnari/documents/Research/FastSummation/fast-bve/include" TYPE FILE FILES
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/amr.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/conservation_fixer.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/direct_sum_utils.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/fast_sum_utils.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/general_utils.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/green_funcs.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/init_utils.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/input_utils.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/interp_utils.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/io_utils.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/mpi_utils.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/rhs_utils.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/structs.hpp"
    "/Users/cygnari/documents/Research/FastSummation/fast-bve/src/vorticity_functions.hpp"
    )
endif()

