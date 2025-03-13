#!/bin/bash
# March 12th 2025, Hongyu Chen

#This script should always be run under the root directory of the project

# Configuration Options
BUILD_TYPE="Release"          # "Debug" or "Release"
ENABLE_MPI="ON"             # "ON" or "OFF"
OPTIMIZATION_LEVEL="-O0"    # Options: "-O0", "-O1", "-O2", "-O3", "-Ofast"
DEBUG_SYMBOLS="OFF"          # "ON" to include debug symbols (-g), "OFF" to exclude

./remove_depend.sh

# Load required modules for NetCDF and HDF5
module load cray-hdf5
module load cray-netcdf

# Define the project root directory (where this script is)
SRC_DIR="$(cd "$(dirname "$0")" && pwd)"

# Remove any in-source CMake artifacts
rm -rf "$SRC_DIR/CMakeCache.txt" "$SRC_DIR/CMakeFiles" "$SRC_DIR/Makefile" "$SRC_DIR/cmake_install.cmake"


if [ -z "$INSTALL_PREFIX" ]; then
  INSTALL_PREFIX="$SRC_DIR"
fi

# Use "./build" as the out-of-source build directory
BUILD_DIR="${SRC_DIR}/build"
rm -rf "$BUILD_DIR"
mkdir "$BUILD_DIR"

DEBUG_FLAGS=""
if [ "$DEBUG_SYMBOLS" == "ON" ]; then
  DEBUG_FLAGS="-g"
fi


cd "$BUILD_DIR" || { echo "Build directory not found: $BUILD_DIR"; exit 1; }

# Configure the project: 
# - The source is pointed to SRC_DIR,
# - The installation prefix is set to SRC_DIR so that install() will copy targets to ${SRC_DIR}/bin.
cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      -DCMAKE_CXX_FLAGS_DEBUG="${OPTIMIZATION_LEVEL} ${DEBUG_FLAGS}" \
      -DENABLE_MPI=${ENABLE_MPI} \
      -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
      "$SRC_DIR"

if [ $? -ne 0 ]; then
  echo "CMake configuration failed. Exiting."
  exit 1
fi

# Build and install the project from the build directory
make && make install

if [ $? -ne 0 ]; then
  echo "Build or installation failed. Exiting."
  exit 1
fi

echo "Build and installation completed successfully."
