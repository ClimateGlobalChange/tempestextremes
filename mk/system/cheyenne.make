# Copyright (c) 2020 Paul Ullrich
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# NCAR Cheyenne

# C++ compiler without and with MPI
CXX=               icpc
MPICXX=            mpicxx

# Additional C++ command line flags
LDFLAGS+= -Wl,-rpath,/glade/u/apps/opt/intel/2017u1/compilers_and_libraries/linux/mkl/lib/intel64

# NetCDF
NETCDF_ROOT=       /glade/u/apps/ch/opt/netcdf/4.4.1.1/intel/17.0.1
NETCDF_CXXFLAGS=   -I$(NETCDF_ROOT)/include
NETCDF_LIBRARIES=  -lnetcdf -lnetcdf_c++
NETCDF_LDFLAGS=    -L$(NETCDF_ROOT)/lib

# DO NOT DELETE
