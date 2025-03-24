# Copyright (c) 2025      Kian Huang, Hongyu Chen
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

list(APPEND SPINE_ARS_FILES
  SpineARs.cpp
)

list(APPEND COMPRESS_BLOBS_FILES
  CompressBlobs.cpp
)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../base
  ${CMAKE_CURRENT_SOURCE_DIR}/../netcdf-cxx-4.2
  ${NetCDF_C_INCLUDE_DIR}
  ${MPI_CXX_INCLUDE_DIRS}
)

add_executable(SpineARs ${SPINE_ARS_FILES})
target_link_libraries(SpineARs PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(CompressBlobs ${COMPRESS_BLOBS_FILES})
target_link_libraries(CompressBlobs PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

install(
  TARGETS SpineARs CompressBlobs
  RUNTIME DESTINATION bin
)
