# Copyright (c) 2025      Kian Huang
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

list(APPEND SPINE_ARS_FILES
  SpineARs.cpp
)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../base
  ${CMAKE_CURRENT_SOURCE_DIR}/../blocking
  ${CMAKE_CURRENT_SOURCE_DIR}/../netcdf-cxx-4.2
  ${NetCDF_C_INCLUDE_DIR}
  ${MPI_CXX_INCLUDE_DIRS}
)

add_executable(SpineARs ${SPINE_ARS_FILES})

target_link_libraries(
  SpineARs PUBLIC
  extremesbase
  netcdf_c++
  NetCDF::NetCDF_C
  ${MPI_CXX_LIBRARIES}
)

install(
  TARGETS SpineARs
)