# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# NERSC Babbage Testbed

CXX=               icpc
F90=               ifort
MPICXX=            mpiicpc
MPIF90=            mpiifort

LDFLAGS+= -Wl,-rpath,/ncar/opt/intel/12.1.0.233/composer_xe_2011_sp1.11.339/mkl/lib/intel64

# NetCDF
NETCDF_ROOT=       /glade/apps/opt/netcdf/4.3.0/intel/12.1.5
NETCDF_CXX_ROOT=   /glade/apps/opt/netcdf/4.3.0/intel/12.1.5
NETCDF_CXXFLAGS=   -I$(NETCDF_ROOT)/include -I$(NETCDF_CXX_ROOT)/include
NETCDF_LIBRARIES=  -lnetcdf -lnetcdf_c++
NETCDF_LDFLAGS=    -L$(NETCDF_ROOT)/lib -L$(NETCDF_CXX_ROOT)/lib -Wl,-rpath=$(NETCDF_CXX_ROOT)/lib

# LAPACK (Intel MKL)
LAPACK_INTERFACE=  FORTRAN
LAPACK_CXXFLAGS=
LAPACK_LIBRARIES=  
LAPACK_LDFLAGS=    -mkl=sequential

# DO NOT DELETE
