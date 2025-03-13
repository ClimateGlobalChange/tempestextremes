# Copyright (c) 2025      Kian Huang, Hongyu Chen
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

list(APPEND ACCUMULATE_DATA_FILES
  AccumulateData.cpp
)

list(APPEND ACCUMULATE_ERA5_FORECAST_FILES
  AccumulateERA5Forecast.cpp
)

list(APPEND AUTO_CURATOR_TOOL_FILES
  AutoCuratorTool.cpp
)

list(APPEND CLIMATOLOGY_FILES
  Climatology.cpp
)

list(APPEND FOURIER_FILTER_FILES
  FourierFilter.cpp
)

list(APPEND GENERATE_CONNECTIVITY_FILES
  GenerateConnectivityFile.cpp
)

list(APPEND GENERATE_NEAREST_NEIGHBOR_MAP_FILES
  GenerateNearestNeighborMap.cpp
)

list(APPEND INTEGRATE_DIMENSION_FILES
  IntegrateDimension.cpp
)

list(APPEND LAGRANGIAN_PARCEL_TRACKER_FILES
  LagrangianParcelTracker.cpp
)

list(APPEND QUANTILE_CALCULATOR_FILES
  QuantileCalculator.cpp
)

list(APPEND SHAPEFILE_MASK_FILES
  ShapefileMask.cpp
)

list(APPEND VARIABLE_PROCESSOR_FILES
  VariableProcessor.cpp
)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../base
  ${CMAKE_CURRENT_SOURCE_DIR}/../blocking
  ${CMAKE_CURRENT_SOURCE_DIR}/../netcdf-cxx-4.2
  ${NetCDF_C_INCLUDE_DIR}
  ${MPI_CXX_INCLUDE_DIRS}
)

add_executable(AccumulateData ${ACCUMULATE_DATA_FILES})
target_link_libraries(AccumulateData PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(AccumulateERA5Forecast ${ACCUMULATE_ERA5_FORECAST_FILES})
target_link_libraries(AccumulateERA5Forecast PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(AutoCuratorTool ${AUTO_CURATOR_TOOL_FILES})
target_link_libraries(AutoCuratorTool PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(Climatology ${CLIMATOLOGY_FILES})
target_link_libraries(Climatology PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(FourierFilter ${FOURIER_FILTER_FILES})
target_link_libraries(FourierFilter PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(GenerateConnectivityFile ${GENERATE_CONNECTIVITY_FILES})
target_link_libraries(GenerateConnectivityFile PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(GenerateNearestNeighborMap ${GENERATE_NEAREST_NEIGHBOR_MAP_FILES})
target_link_libraries(GenerateNearestNeighborMap PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(IntegrateDimension ${INTEGRATE_DIMENSION_FILES})
target_link_libraries(IntegrateDimension PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(LagrangianParcelTracker ${LAGRANGIAN_PARCEL_TRACKER_FILES})
target_link_libraries(LagrangianParcelTracker PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(QuantileCalculator ${QUANTILE_CALCULATOR_FILES})
target_link_libraries(QuantileCalculator PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(ShapefileMask ${SHAPEFILE_MASK_FILES})
target_link_libraries(ShapefileMask PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(VariableProcessor ${VARIABLE_PROCESSOR_FILES})
target_link_libraries(VariableProcessor PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

install(
  TARGETS
    AccumulateData
    AccumulateERA5Forecast
    AutoCuratorTool
    Climatology
    FourierFilter
    GenerateConnectivityFile
    GenerateNearestNeighborMap
    IntegrateDimension
    LagrangianParcelTracker
    QuantileCalculator
    ShapefileMask
    VariableProcessor
  RUNTIME DESTINATION bin
)

