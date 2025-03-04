# Copyright (c) 2025      Kian Huang
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.12)

project(tempestextremes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE ${CMAKE_CURRENT_SOURCE_DIR})
endif()

option(CMAKE_BUILD_TYPE "Build type" Release)
option(ENABLE_MPI "Enable MPI" ON)

if(WIN32)
  add_compile_definitions(_USE_MATH_DEFINES 1)
endif()

if(ENABLE_MPI)
  find_package(MPI)

  if(MPI_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTEMPEST_MPIOMP")
  endif()
endif()

# Required dependencies
find_package(NetCDF REQUIRED)

add_subdirectory(src)