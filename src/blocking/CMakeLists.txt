# Copyright (c) 2025      Kian Huang, Hongyu Chen
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

list(APPEND BLOCKING_UTILITIES_FILES
  BlockingUtilities.cpp
)

list(APPEND AVG_VAR_FILES
  ${BLOCKING_UTILITIES_FILES}
  AvgVar.cpp
)

list(APPEND BLOCKING_DEVS_FILES
  ${BLOCKING_UTILITIES_FILES}
  BlockingDevs.cpp
)

list(APPEND BLOCKING_DFT_FILES
  ${BLOCKING_UTILITIES_FILES}
  BlockingDFT.cpp
  DFT.cpp
)

list(APPEND BLOCKING_GHG_FILES
  ${BLOCKING_UTILITIES_FILES}
  BlockingGHG.cpp
)

list(APPEND BLOCKING_NORM_DEVS_FILES
  ${BLOCKING_UTILITIES_FILES}
  BlockingNormDevs.cpp
)

list(APPEND BLOCKING_PV_FILES
  ${BLOCKING_UTILITIES_FILES}
  BlockingPV.cpp
  Interpolate.cpp
)

list(APPEND BLOCKING_THRESH_FILES
  ${BLOCKING_UTILITIES_FILES}
  BlockingThresh.cpp
  DFT.cpp
)

list(APPEND COMBINE_BLOBS_FILES
  ${BLOCKING_UTILITIES_FILES}
  CombineBlobs.cpp
)

list(APPEND DIALY_AVERAGE_FILES
  ${BLOCKING_UTILITIES_FILES}
  DailyAverage.cpp
)

list(APPEND DENSITY_CALCULATIONS_FILES
  ${BLOCKING_UTILITIES_FILES}
  DensityCalculations.cpp
)

list(APPEND DETREND_HEIGHTS_FILES
  ${BLOCKING_UTILITIES_FILES}
  DetrendHeights.cpp
)

list(APPEND EXTRACT_TIME_STEP_FILES
  ${BLOCKING_UTILITIES_FILES}
  ExtractTimeStep.cpp
)

list(APPEND SMOOTH_61_DAY_FILES
  ${BLOCKING_UTILITIES_FILES}
  Smooth61Day.cpp
)

list(APPEND SPLIT_FILE_FILES
  ${BLOCKING_UTILITIES_FILES}
  SplitFile.cpp
)

list(APPEND VAR4D_TO_3D_FILES
  ${BLOCKING_UTILITIES_FILES}
  Interpolate.cpp
  Var4Dto3D.cpp
)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../base
  ${CMAKE_CURRENT_SOURCE_DIR}/../blocking
  ${CMAKE_CURRENT_SOURCE_DIR}/../netcdf-cxx-4.2
  ${NetCDF_C_INCLUDE_DIR}
  ${MPI_CXX_INCLUDE_DIRS}
)

add_executable(AvgVar ${AVG_VAR_FILES})
target_link_libraries(AvgVar PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(BlockingDevs ${BLOCKING_DEVS_FILES})
target_link_libraries(BlockingDevs PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(BlockingDFT ${BLOCKING_DFT_FILES})
target_link_libraries(BlockingDFT PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(BlockingGHG ${BLOCKING_GHG_FILES})
target_link_libraries(BlockingGHG PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(BlockingNormDevs ${BLOCKING_NORM_DEVS_FILES})
target_link_libraries(BlockingNormDevs PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(BlockingPV ${BLOCKING_PV_FILES})
target_link_libraries(BlockingPV PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(BlockingThresh ${BLOCKING_THRESH_FILES})
target_link_libraries(BlockingThresh PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(CombineBlobs ${COMBINE_BLOBS_FILES})
target_link_libraries(CombineBlobs PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(DailyAverage ${DIALY_AVERAGE_FILES})
target_link_libraries(DailyAverage PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(DensityCalculations ${DENSITY_CALCULATIONS_FILES})
target_link_libraries(DensityCalculations PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(DetrendHeights ${DETREND_HEIGHTS_FILES})
target_link_libraries(DetrendHeights PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(ExtractTimeStep ${EXTRACT_TIME_STEP_FILES})
target_link_libraries(ExtractTimeStep PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(Smooth61Day ${SMOOTH_61_DAY_FILES})
target_link_libraries(Smooth61Day PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(SplitFile ${SPLIT_FILE_FILES})
target_link_libraries(SplitFile PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})

add_executable(Var4Dto3D ${VAR4D_TO_3D_FILES})
target_link_libraries(Var4Dto3D PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})


install(
  TARGETS
    AvgVar
    BlockingDevs
    BlockingDFT
    BlockingGHG
    BlockingNormDevs
    BlockingPV
    BlockingThresh
    CombineBlobs
    DailyAverage
    DensityCalculations
    DetrendHeights
    ExtractTimeStep
    Smooth61Day
    SplitFile
    Var4Dto3D
  RUNTIME DESTINATION bin
)