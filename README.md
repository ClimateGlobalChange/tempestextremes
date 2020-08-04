TempestExtremes
================

Author:  Paul Ullrich
Email:   paullrich@ucdavis.edu

Copyright 2020 Paul Ullrich

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Installation
============
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

Usage
=====
Details of the various executables that are part of TempestExtremes can be found in the user guide:
https://climate.ucdavis.edu/tempestextremes.php

Publications
============
If you use the TempestExtremes software please cite our publications:

[http://dx.doi.org/10.5194/gmd-2016-217] Ullrich, P.A. and C.M. Zarzycki (2017) "TempestExtremes v1.0: A framework for scale-insensitive pointwise feature tracking on unstructured grids" Geosci. Model. Dev. 10, pp. 1069-1090, doi: 10.5194/gmd-10-1069-2017. 

[http://dx.doi.org/10.1002/2016GL071606] Zarzycki, C.M. and P.A. Ullrich (2017) "Assessing sensitivities in algorithmic detection of tropical cyclones in climate data" Geophys. Res. Lett. 44 (2), pp. 1141-1149, doi: 10.1002/2016GL071606. 
