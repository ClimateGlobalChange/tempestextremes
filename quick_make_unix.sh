#!/bin/bash
# 2025, Hongyu Chen




# This script should always be run under the root directory of the project.
# It provides a general, quick command to build the project.
# Please ensure the required NetCDF and HDF5 (and MPI, if needed) are available.

# Configuration Options
BUILD_TYPE="Release"          # "Debug" or "Release"
ENABLE_MPI="ON"               # "ON" or "OFF"
OPTIMIZATION_LEVEL="-O0"      # Options: "-O0", "-O1", "-O2", "-O3", "-Ofast"
DEBUG_SYMBOLS="OFF"           # "ON" to include debug symbols (-g), "OFF" to exclude
INSTALL_PREFIX=""             # Specify the installation directory.
                              # If left blank, it defaults to the project root (TEMPEST_EXTREMES_SOURCE_DIR)
                              # and final executables will be installed in TEMPEST_EXTREMES_SOURCE_DIR/bin.

./remove_depend.sh



# Detect system type
OS_TYPE="$(uname -s)"
SYSTEM_TYPE=""

if [ "$OS_TYPE" == "Darwin" ]; then
    SYSTEM_TYPE="MacOS/Linux"
elif [ "$OS_TYPE" == "Linux" ]; then
    # Check for NERSC Perlmutter: expect NERSC_HOST set to "perlmutter"
    if [ -n "$NERSC_HOST" ] && [ "$NERSC_HOST" = "perlmutter" ]; then
        SYSTEM_TYPE="NERSC Perlmutter"
    # Check if hostname contains "derecho" (case-insensitive)
    elif echo "$HOSTNAME" | grep -qi "derecho"; then
        SYSTEM_TYPE="NCAR Derecho"
    else
        SYSTEM_TYPE="MacOS/Linux"
    fi
elif [[ "$OS_TYPE" == *"CYGWIN"* || "$OS_TYPE" == *"MINGW"* || "$OS_TYPE" == *"MSYS"* ]]; then
    SYSTEM_TYPE="Windows"
else
    SYSTEM_TYPE="Unknown"
fi

echo "Detected system: ${SYSTEM_TYPE}"

# System-specific module loading
if [ "$SYSTEM_TYPE" = "MacOS/Linux" ]; then
    echo "Running on MacOS/Linux. Proceeding with default configuration..."
    # No additional module load commands needed.
elif [ "$SYSTEM_TYPE" = "NERSC Perlmutter" ]; then
    echo "Loading modules for NERSC Perlmutter..."
    module load cray-hdf5
    module load cray-netcdf
elif [ "$SYSTEM_TYPE" = "NCAR Derecho" ]; then
    echo "Loading modules for NCAR Derecho..."
    module load cmake
    module load ncarenv
    module load ncarcompilers
    module load intel
    module load cray-mpich
    module load netcdf
elif [ "$SYSTEM_TYPE" = "Windows" ]; then
    echo "Windows detected. Please follow the README instructions for Windows build or manually run the commands in your bash enviroment."
    exit 1
else
    echo "Unable to detect the system. Please refer to the README instructions or mannually try commands in ./quick_make_general.sh."
    exit 1
fi

./remove_depend.sh

# Load required modules for NetCDF and HDF5
module load cray-hdf5
module load cray-netcdf

# Define the project root directory (where this script is)
SRC_DIR="$(cd "$(dirname "$0")" && pwd)"

# Remove any in-source CMake artifacts
rm -rf "$SRC_DIR/CMakeCache.txt" "$SRC_DIR/CMakeFiles" "$SRC_DIR/Makefile" "$SRC_DIR/cmake_install.cmake"

# Set INSTALL_PREFIX to SRC_DIR if not provided
if [ -z "$INSTALL_PREFIX" ]; then
  INSTALL_PREFIX="$SRC_DIR"
fi

# Clean up the installed binary directory (./bin) before building
rm -rf "$INSTALL_PREFIX/bin"

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
# - The source directory is set to SRC_DIR.
# - The installation prefix is set to INSTALL_PREFIX so that install() will copy targets to ${INSTALL_PREFIX}/bin.
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

# (Optional) Cleanup step:
# The build directory (./build) is used for development and debugging.
# It contains all intermediate build files and temporary artifacts.
# The final, user-deliverable executables are installed to ./bin.
# It is not recommended to mix these directories.

# For end users who want a clean structure, you can remove the build directory.
# Developers or those debugging might prefer to keep it for faster incremental builds.
make clean
echo "Cleaned up the ${SRC_DIR}/build directory. All executables are located in ${INSTALL_PREFIX}/bin."

