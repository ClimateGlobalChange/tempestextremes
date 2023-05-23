# Copyright (c) 2023 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# NERSC Perlmutter

# C++ compiler without and with MPI
CXX=               CC
MPICXX=            CC

# NetCDF C library arguments
NETCDF_ROOT=       $(NETCDF_DIR)
NETCDF_CXXFLAGS=   -I$(NETCDF_ROOT)/include
NETCDF_LIBRARIES=  -lnetcdf -lnetcdf_c++
NETCDF_LDFLAGS=    -L$(NETCDF_ROOT)/lib

# DO NOT DELETE
