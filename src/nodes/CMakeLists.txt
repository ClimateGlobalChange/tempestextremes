# Copyright (c) 2025      Kian Huang, Hongyu Chen
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

list(APPEND DETECT_NODES_FILES
  DetectNodes.cpp
)

list(APPEND HISTOGRAM_NODES_FILES
  HistogramNodes.cpp
)

list(APPEND NODE_FILE_COMPOSE_FILES
  NodeFileCompose.cpp
)

list(APPEND NODE_FILE_EDITOR_FILES
  NodeFileEditor.cpp
)

list(APPEND NODE_FILE_FILTER_FILES
  NodeFileFilter.cpp
)

list(APPEND STITCH_NODES_FILES
  StitchNodes.cpp
)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../base
  ${CMAKE_CURRENT_SOURCE_DIR}/../netcdf-cxx-4.2
  ${NetCDF_C_INCLUDE_DIR}
  ${MPI_CXX_INCLUDE_DIRS}
)

add_executable(DetectNodes ${DETECT_NODES_FILES})
target_link_libraries(DetectNodes PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(HistogramNodes ${HISTOGRAM_NODES_FILES})
target_link_libraries(HistogramNodes PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(NodeFileCompose ${NODE_FILE_COMPOSE_FILES})
target_link_libraries(NodeFileCompose PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(NodeFileEditor ${NODE_FILE_EDITOR_FILES})
target_link_libraries(NodeFileEditor PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(NodeFileFilter ${NODE_FILE_FILTER_FILES})
target_link_libraries(NodeFileFilter PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(StitchNodes ${STITCH_NODES_FILES})
target_link_libraries(StitchNodes PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

install(
  TARGETS DetectNodes HistogramNodes NodeFileCompose NodeFileEditor NodeFileFilter StitchNodes
  RUNTIME DESTINATION bin
)
