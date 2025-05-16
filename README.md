# TempestExtremes


Author:  Paul Ullrich
Email:   paullrich@ucdavis.edu

Copyright 2025 Paul Ullrich

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Installation via conda

TempestExtremes can be found on conda-forge here:

https://anaconda.org/conda-forge/tempest-extremes

To install from conda use the command line:

```conda install -c conda-forge tempest-extremes```


# Installation via CMake (Recommended)

TempestExtremes can be built and installed on various systems using CMake. Our new script, `./quick_make_unix.sh` in the TempestExtremes root directory, automatically detects your platform(UNIX-based systems only) and loads any required modules before building. To install all TempestExtremes executable to the ./bin directory, just run the quick make script.

## General CMake Configuration

- **Install Prefix:** Set with `-DCMAKE_INSTALL_PREFIX=PATH_TO_INSTALL`. By default, it installs to the project root, with executables placed in `TEMPEST_EXTREMES_SOURCE_DIR/bin`.
- **Build Type:** Specify via `-DCMAKE_BUILD_TYPE=[Release/Debug]`.
- **MPI Support:** Enable or disable using `-DENABLE_MPI=ON` or `-DENABLE_MPI=OFF`.
- **Out-of-Source Build:** Build files are generated in `./build/bin` (keeping the source directory clean).
- **Installation Locations:** Final executables are installed to `./bin` and libraries/archives to `./lib`.

## Quick-Make Scripts for generic Linux/Unix based systems
A ready-to-use `./quick_make_unix.sh` script is available for end users to run on common UNIX-based platforms such as MacOS and Linux, as well as on command HPC systems like NERSC Perlmutter and NCAR Derecho.
### System/Platforms Detection in `quick_make_unix.sh`
- **MacOS/Linux (Generic):** The script runs the default commands.
- **NERSC Perlmutter:** The script automatically loads `cray-hdf5` and `cray-netcdf` modules before building.
- **NCAR Derecho:** The script loads required modules such as `ncarenv`, `ncarcompilers`, `intel`, `cray-mpich`, and `netcdf`.
- **Windows:** This bash script is prepared for UNIX-like environments. If you're on Windows, please refer to the Windows instructions below. If you have a bash environment on Windows but the script fails to run automatically, simply run the commands manually starting from `./remove_depend.sh`.
- **Unknown/Other:** If your system isn’t recognized or errors occur (often due to non-standard dependency paths), don’t panic—simply run the commands manually starting from `./remove_depend.sh` and fix the related errors shown in the terminal.


Notes:
- **For End Users:**  
  If you only need the final deliverables, run the quick-make script. It will install executables to `./bin` and (by default) remove the build directory for a clean structure.

- **For Developers:**  
  The `./build` directory is used for development and debugging. It contains all intermediate files, which helps speed up incremental builds. Developers should keep the build directory intact (by commenting out the cleanup step) for faster rebuilds.  


To use `./quick_make_unix.sh` , update the configuration options in the script:
```bash
# Configuration Options
BUILD_TYPE="Release"          # "Debug" or "Release"
ENABLE_MPI="ON"               # "ON" or "OFF"
OPTIMIZATION_LEVEL="-O3"      # Options: "-O0", "-O1", "-O2", "-O3", "-Ofast"
DEBUG_SYMBOLS="OFF"           # "ON" to include debug symbols (-g), "OFF" to exclude
INSTALL_PREFIX=""             # Specify the installation directory. If left blank, it defaults to 
                              # the project root (TEMPEST_EXTREMES_SOURCE_DIR) and final executables 
                              # will be installed in TEMPEST_EXTREMES_SOURCE_DIR/bin.
```
The run `./quick_make_unix.sh`.

## Unix/Linux-Based Systems (Manual Install)

If you are using Nersc Permultter, you will need to run this first:
```
module load cray-hdf5
module load cray-netcdf
```
If you are using NCAR Derecho, you will need to run this first:
```
module load cmake ncarenv gcc craype cray-mpich hdf5 netcdf
```
Use the following commands to compile on Unix- or Linux-based systems ):
```
cd TEMPEST_EXTREMES_SOURCE_DIR
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=[Release/Debug] -DCMAKE_INSTALL_PREFIX=PATH_TO_INSTALL -DCMAKE_PREFIX_PATH=MPI_ROOT ..
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


## Developer Notes


1. **Using the Quick Make Script:**  
   If you are a developer using `./quick_make_general.sh`, remember to comment out or remove the cleanup step at the end (e.g. `make clean`). It's not efficient to rebuild the entire project every time you make changes—this script is mainly for a full rebuild or for end users who want a clean install. For faster incremental builds during development, build manually so that intermediate files are preserved.

2. **Adding New Executables:**  
   To add a new executable to the project, refer to following examples:

   **Example: Adding a New Executable**  
   For example, to add `NewBlobsFeature`:
   1. **Edit the CMakeLists.txt File:**  
      Open the appropriate file (e.g., `src/blobs/CMakeLists.txt`) and add:
      ```cmake
      add_executable(NewBlobsFeature NewBlobsFeature.cpp)
      target_link_libraries(NewBlobsFeature PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})
      ```
   2. **Update the Install Rule:**  
      Modify the install command in the same CMakeLists.txt to include the new executable:
      ```cmake
      install(
        TARGETS BlobStats DetectBlobs PersistentBlobs StitchBlobs NewBlobsFeature
        RUNTIME DESTINATION bin
      )
      ```
      This ensures that `NewBlobsFeature` will be copied to `./bin` when you run `make install`.

3. **Combining Multiple Files into One Executable:**  
   If your new executable requires multiple source files, do the following:

   **Example: Combining Multiple Files into One Executable**  
   For instance, to create an executable named `BlockingUtilities` from multiple files:
   1. **Edit the CMakeLists.txt File:**  
      Open `src/blocking/CMakeLists.txt` and add:
      ```cmake
      list(APPEND BLOCKINGUTILITIES
        BlockingUtilities.h
        BlockingUtilities.cpp
        MoreBlockingUtilities.cpp
      )

      add_executable(BlockingUtilities ${BLOCKINGUTILITIES})
      target_link_libraries(BlockingUtilities PUBLIC extremesbase netcdf_c++ ${MPI_CXX_LIBRARIES})
      ```
   2. **Update the Install Rule:**  
      Ensure that the new target is included in the install command:
      ```cmake
      install(
        TARGETS BlockingUtilities
        RUNTIME DESTINATION bin
      )
      ```

By following these guidelines, you can efficiently develop and extend TempestExtremes while keeping a clean separation between the build artifacts and the final deliverables.



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
