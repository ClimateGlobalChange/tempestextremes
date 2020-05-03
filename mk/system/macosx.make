# Copyright (c) 2020 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Mac OSX (Darwin)

# C++ compiler without and with MPI
CXX=               g++
MPICXX=            mpic++

# Additional C++ command line flags

# NetCDF C library arguments
NETCDF_ROOT=       /opt/local
NETCDF_CXXFLAGS=   -I$(NETCDF_ROOT)/include
NETCDF_LIBRARIES=  -lnetcdf -lnetcdf_c++
NETCDF_LDFLAGS=    -L$(NETCDF_ROOT)/lib -Wl,-rpath,$(NETCDF_CXX_ROOT)/lib

# DO NOT DELETE
