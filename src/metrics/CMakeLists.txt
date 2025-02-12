# Copyright (c) 2025      Kian Huang
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

list(APPEND HEAT_WAVE_METRICS_FILES
  HeatWaveMetrics.cpp
)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../base
  ${CMAKE_CURRENT_SOURCE_DIR}/../blocking
  ${CMAKE_CURRENT_SOURCE_DIR}/../netcdf-cxx-4.2
  ${NetCDF_C_INCLUDE_DIR}
  ${MPI_CXX_INCLUDE_DIRS}
)

add_executable(HeatWaveMetrics ${HEAT_WAVE_METRICS_FILES})
target_link_libraries(HeatWaveMetrics PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

install(
  TARGETS HeatWaveMetrics
)