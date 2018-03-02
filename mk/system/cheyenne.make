# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# NCAR Cheyenne

CXX=               icpc
F90=               ifort
MPICXX=            mpicxx
MPIF90=            mpif90

LDFLAGS+= -Wl,-rpath,/glade/u/apps/opt/intel/2017u1/compilers_and_libraries/linux/mkl/lib/intel64

# NetCDF
NETCDF_ROOT=       /glade/u/apps/ch/opt/netcdf/4.4.1.1/intel/17.0.1
NETCDF_CXX_ROOT=   /glade/u/apps/ch/opt/netcdf/4.4.1.1/intel/17.0.1
NETCDF_CXXFLAGS=   -I$(NETCDF_ROOT)/include -I$(NETCDF_CXX_ROOT)/include
NETCDF_LIBRARIES=  -lnetcdf -lnetcdf_c++
NETCDF_LDFLAGS=    -L$(NETCDF_ROOT)/lib -L$(NETCDF_CXX_ROOT)/lib -Wl,-rpath=$(NETCDF_CXX_ROOT)/lib

# LAPACK (Intel MKL)
LAPACK_INTERFACE=  FORTRAN
LAPACK_CXXFLAGS=
LAPACK_LIBRARIES=  
LAPACK_LDFLAGS=    -mkl=sequential

# DO NOT DELETE
