# Copyright (c) 2025      Kian Huang, Hongyu Chen
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

list(APPEND FILES
  Announce.h
  Announce.cpp
  ArgumentTree.h
  ArgumentTree.cpp
  AutoCurator.h
  AutoCurator.cpp
  CommandLine.h
  Constants.h
  CoordTransforms.h
  DataArray1D.h
  DataArray2D.h
  DataArray3D.h
  DataArray4D.h
  DataOp.h
  DataOp.cpp
  Defines.h
  Exception.h
  FilenameList.h
  FiniteElementTools.h
  FiniteElementTools.cpp
  FourierTransforms.h
  FunctionTimer.h
  FunctionTimer.cpp
  GaussLobattoQuadrature.h
  GaussLobattoQuadrature.cpp
  GaussQuadrature.h
  GaussQuadrature.cpp
  GridElements.h
  GridElements.cpp
  kdtree.h
  kdtree.cpp
  LatLonBox.h
  LegendrePolynomial.h
  LegendrePolynomial.cpp
  lodepng.h
  lodepng.cpp
  MathExpression.h
  MeshUtilities.h
  MeshUtilities.cpp
  MeshUtilitiesFuzzy.h
  MeshUtilitiesFuzzy.cpp
  NcFileVector.h
  NcFileVector.cpp
  NetCDFUtilities.h
  NetCDFUtilities.cpp
  NodeFileUtilities.h
  NodeFileUtilities.cpp
  order32.h
  PolynomialInterp.h
  PolynomialInterp.cpp
  RLLPolygonArray.h
  RLLPolygonArray.cpp
  ShpFile.h
  ShpFile.cpp
  SimpleGrid.h
  SimpleGrid.cpp
  SimpleGridUtilities.h
  SimpleGridUtilities.cpp
  SparseMatrix.h
  STLStringHelper.h
  Subscript.h
  ThresholdOp.h
  ThresholdOp.cpp
  TimeMatch.h
  TimeObj.h
  TimeObj.cpp
  Units.h
  Variable.h
  Variable.cpp
)

add_library(
  extremesbase
  STATIC
  ${FILES}
)

target_include_directories(extremesbase
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../netcdf-cxx-4.2
  ${NetCDF_C_INCLUDE_DIR}
  ${MPI_CXX_INCLUDE_DIRS}
)

set_target_properties(
  extremesbase PROPERTIES
  LINKER_LANGUAGE CXX
)

# Install the static library to the "./lib" folder
install(
  TARGETS extremesbase
  ARCHIVE DESTINATION lib
)
