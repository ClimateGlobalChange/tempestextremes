# Copyright (c) 2020 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# UC Davis Agri 

# C++ compiler without and with MPI
CXX=               g++
MPICXX=            mpiCC

# Additional C++ command line flags
CXXFLAGS+=         -fPIC

# NetCDF C library arguments
NETCDF_ROOT=       $(NETCDF_HOME)
NETCDF_CXXFLAGS=   -I$(NETCDF_ROOT)/include
NETCDF_LIBRARIES=  -lnetcdf_c++ -lnetcdf
NETCDF_LDFLAGS=    -L$(NETCDF_ROOT)/lib

# DO NOT DELETE
