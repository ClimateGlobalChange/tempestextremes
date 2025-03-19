# Copyright (c) 2025      Kian Huang, Hongyu Chen
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.12)

project(tempestextremes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

# Do not override the installation prefix to the source directory
# Users should specify an install prefix if desired via -DCMAKE_INSTALL_PREFIX

# Set default build type if not already defined
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# MPI is now enabled only via the command line (-DENABLE_MPI=ON)
if(DEFINED ENABLE_MPI AND ENABLE_MPI)
  find_package(MPI REQUIRED)
  if(MPI_FOUND)
    add_compile_definitions(TEMPEST_MPIOMP)
  endif()
endif()

if(WIN32)
  add_compile_definitions(_USE_MATH_DEFINES 1)
endif()

# Required dependencies
find_package(NetCDF REQUIRED)

# Output directories for out-of-source builds
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

add_subdirectory(src)

# Ensure the runtime output directory exists.
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Install rule: copy all files from the build's bin directory to the install prefix's bin directory.
install(DIRECTORY "${CMAKE_BINARY_DIR}/bin/" DESTINATION bin)