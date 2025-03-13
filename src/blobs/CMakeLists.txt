# Copyright (c) 2025      Kian Huang, Hongyu Chen
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

list(APPEND BLOB_STATS_FILES
  BlobStats.cpp
)

list(APPEND DETECT_BLOBS_FILES
  DetectBlobs.cpp
)

list(APPEND PERSISTENT_BLOBS_FILES
  PersistentBlobs.cpp
)

list(APPEND STITCH_BLOBS_FILES
  StitchBlobs.cpp
)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../base
  ${CMAKE_CURRENT_SOURCE_DIR}/../blocking
  ${CMAKE_CURRENT_SOURCE_DIR}/../netcdf-cxx-4.2
  ${NetCDF_C_INCLUDE_DIR}
  ${MPI_CXX_INCLUDE_DIRS}
)

add_executable(BlobStats ${BLOB_STATS_FILES})
target_link_libraries(BlobStats PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(DetectBlobs ${DETECT_BLOBS_FILES})
target_link_libraries(DetectBlobs PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(PersistentBlobs ${PERSISTENT_BLOBS_FILES})
target_link_libraries(PersistentBlobs PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(StitchBlobs ${STITCH_BLOBS_FILES})
target_link_libraries(StitchBlobs PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

# Install executables to the "bin" directory under the installation prefix.
install(
  TARGETS BlobStats DetectBlobs PersistentBlobs StitchBlobs
  RUNTIME DESTINATION bin
)