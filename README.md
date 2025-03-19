TempestExtremes
================

Author:  Paul Ullrich
Email:   paullrich@ucdavis.edu

Copyright 2025 Paul Ullrich

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Installation via conda
======================
TempestExtremes can be found on conda-forge here:

https://anaconda.org/conda-forge/tempest-extremes

To install from conda use the command line:

```conda install -c conda-forge tempest-extremes```

Installation via make
=====================
TempestExtremes uses a basic make-based build system that has been tested in multiple Unix and Linux based environments.  TempestExtremes currently only requires the NetCDF C library provided as an external dependency.  Linker commands needed to include -lnetcdf must be specified in a system makefile found in mk/system/.  Once the compiler executable and flags have been set, executing "make" in the tempestextremes directory will perform the build, storing binaries in the "bin" directory.

To compile on NERSC Cori:
Execute "module load cray-netcdf" prior to compilation.  Compiler executable and flags are specified in mk/system/cori.make.

To compile on NCAR systems, including Cheyenne and Casper:
No modifications are needed.  Compiler executable and flags are specified in mk/system/cheyenne.make.

To prepare compilation on MacOSX:
Modify mk/system/macosx.make to specify desired compiler executable and flags.

To prepare compilation on another Unix- or Linux-based system:
Modify mk/system/default.make to specify appropriate compiler executable and flags.

Parallel compilation with MPI is enabled by default.  To compile TempestExtremes as a serial product edit "mk/config.mk" and change "PARALLEL= MPIOMP" to "PARALLEL=NONE".

Note:  If the build fails or is cancelled part way through the process, it may be necessary to remove extraneous dependency files.  To do so run "remove_depend.sh" prior to recompiling.

Installation via CMake
=====================
Using CMake, one can easily install TempestExtremes on different operational systems with the required compiler and dependencies. Dependencies can be downloaded and installed manually or via package management software (e.g., [Conda](https://docs.conda.io)) and environment variables can be set to help CMake locate those dependencies (see CMake/FindNetCDF.cmake and [CMake help](https://cmake.org/cmake/help/latest/module/FindMPI.html)).

General CMake configuration details:
- **Install Prefix:** Specify an installation prefix via `-DCMAKE_INSTALL_PREFIX=PATH_TO_INSTALL` if desired. The default will be at the source directory of the project.
- **Build Type:** Manually set the build type ("Release" or "Debug") via `-DCMAKE_BUILD_TYPE=[Release/Debug]`.
- **MPI Enable:** Manually enable or disable MPI support using `-DENABLE_MPI=ON` or `-DENABLE_MPI=OFF`.
- **Out-of-Source Build:** For best practices, build files are written to `./build/bin` by default.
- **Installation Locations:** Executables are installed to `./bin` and libraries/archives to `./lib`.

Notes:
- **For End Users:**  
  If you are an end user and want a clean structure with only the final deliverables, simply run the provided "quick make" scripts. The executables will be copied to `./bin`, and the script will remove the build directory after installation by default.

- **For Developers:**  
  The `./build` directory is used for development and debugging—it holds all intermediate build files, which helps with faster incremental builds. It is recommended to keep the build directory intact during development and comment out any cleanup steps that remove it if you are running one of those provided "quick make" scripts.

Quick-Make Scripts
------------------
Two scripts are provided:
- `./quick_make_general.sh` for general systems.
- `./quick_make_perlmutter.sh` for NERSC Perlmutter.

To use them, update the configuration options in the script:
```bash
# Configuration Options
BUILD_TYPE="Release"          # "Debug" or "Release"
ENABLE_MPI="ON"               # "ON" or "OFF"
OPTIMIZATION_LEVEL="-O0"      # Options: "-O0", "-O1", "-O2", "-O3", "-Ofast"
DEBUG_SYMBOLS="OFF"           # "ON" to include debug symbols (-g), "OFF" to exclude
INSTALL_PREFIX=""             # Specify the installation directory. If left blank, it defaults to 
                              # the project root (TEMPEST_EXTREMES_SOURCE_DIR) and final executables 
                              # will be installed in TEMPEST_EXTREMES_SOURCE_DIR/bin.
```
The run `./quick_make_general.sh` or `./quick_make_perlmutter.sh` for NERSC Perlmutter.

## Unix/Linux-Based Systems
Use the following commands to compile on Unix- or Linux-based systems ([netCDF](https://downloads.unidata.ucar.edu/netcdf/) required):
```
cd TEMPEST_EXTREMES_SOURCE_DIR
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=[Release/Debug] -DCMAKE_INSTALL_PREFIX=PATH_TO_INSTALL -DNetCDF_PATH=NETCDF_ROOT -DCMAKE_PREFIX_PATH=MPI_ROOT ..
make && make install
```

## Windows (Experimental)
Follow these steps to compile on the Windows system ([Visual Studio](https://visualstudio.microsoft.com) and [CMake](https://cmake.org/download/) required):
1. Download and install [netCDF](https://downloads.unidata.ucar.edu/netcdf/) and [Microsoft MPI](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi). Select "Add netCDF to the system PATH for the current user" when instaling netCDF.
2. Use the following PowerShell commands to generate a Visual Studio project.
```
$Env:MSMPI_LIB64="MSMPI_SDK_ROOT/Lib/x64"
$Env:MSMPI_INC="MSMPI_SDK_ROOT/Include"
cd TEMPEST_EXTREMES_SOURCE_DIR
mkdir build
cd build
cmake -G "Visual Studio 17" -DCMAKE_INSTALL_PREFIX=PATH_TO_INSTALL -DNetCDF_ROOT=YOUR_NETCDF_INSTALLATION_PATH ..
```
3. Open the `tempestextremes.sln` file in the `./build` directory with Visual Studio and Build the software by clicking "Build" and "Build INSTALL" at the top menu bar. [Learn more about building Visual Studio projects](https://learn.microsoft.com/en-us/visualstudio/ide/building-and-cleaning-projects-and-solutions-in-visual-studio). 

## HPC Systems

### NERSC Perlmutter
Use the following commands to compile on [Perlmutter](https://docs.nersc.gov/systems/perlmutter/running-jobs/):
```
module load cray-hdf5
module load cray-netcdf

cd TEMPEST_EXTREMES_SOURCE_DIR
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=[Release/Debug] -DCMAKE_INSTALL_PREFIX=PATH_TO_INSTALL ..
make && make install
```
Additionally, a ready-to-run script (`./quick_make_perlmutter.sh`) is provided for Perlmutter. Simply update its configuration options and run `./quick_make_perlmutter.sh`

### NCAR Derecho
Use the following commands to compile on [Derecho](https://ncar-hpc-docs.readthedocs.io/en/latest/compute-systems/derecho/compiling-code-on-derecho/):
```
module load ncarenv/23.09
module load ncarcompilers/1.0.0
module load intel/2023.2.1
module load cray-mpich/8.1.27
module load netcdf/4.9.2
cd TEMPEST_EXTREMES_SOURCE_DIR
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=[Release/Debug] -DCMAKE_INSTALL_PREFIX=PATH_TO_INSTALL ..
make && make install
```

Usage
=====
Details of the various executables that are part of TempestExtremes can be found in the user guide:
https://climate.ucdavis.edu/tempestextremes.php

Publications
============
If you use the TempestExtremes software please cite our publications:

[https://dx.doi.org/10.5194/gmd-14-5023-2021] Ullrich, P.A., C.M. Zarzycki, E.E. McClenny, M.C. Pinheiro, A.M. Stansfield and K.A. Reed (2021) "TempestExtremes v2.1: A community framework for feature detection, tracking and analysis in large datasets" Geosci. Model. Dev. 14, pp. 5023–5048, doi: 10.5194/gmd-14-5023-2021.

[http://dx.doi.org/10.5194/gmd-2016-217] Ullrich, P.A. and C.M. Zarzycki (2017) "TempestExtremes v1.0: A framework for scale-insensitive pointwise feature tracking on unstructured grids" Geosci. Model. Dev. 10, pp. 1069-1090, doi: 10.5194/gmd-10-1069-2017. 

[http://dx.doi.org/10.1002/2016GL071606] Zarzycki, C.M. and P.A. Ullrich (2017) "Assessing sensitivities in algorithmic detection of tropical cyclones in climate data" Geophys. Res. Lett. 44 (2), pp. 1141-1149, doi: 10.1002/2016GL071606. 
