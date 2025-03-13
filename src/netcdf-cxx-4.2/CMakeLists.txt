# Copyright (c) 2025      Kian Huang, Hongyu Chen
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

list(APPEND FILES
  ncvalues.cpp
  netcdf.cpp
)

add_library(
  netcdf_c++
  STATIC
  ${FILES}
)

target_include_directories(
  netcdf_c++
  PRIVATE
  ${NetCDF_C_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(
  netcdf_c++
  PUBLIC
  NetCDF::NetCDF_C
  ${MPI_CXX_LIBRARIES}
)

# Install the static library to the "lib" folder
install(
  TARGETS netcdf_c++
  ARCHIVE DESTINATION lib
)
