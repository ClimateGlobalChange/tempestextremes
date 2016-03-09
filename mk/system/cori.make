# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

CXX= 		g++
F90=		gfortran
MPICXX=		mpiCC
MPIF90=		mpif90

CXXFLAGS+=	-fPIC
F90FLAGS+=	-fPIC
#LIBRARIES+=
#LDFLAGS+=	

F90_RUNTIME=	-lgfortran

# HPX
#HPX_CXXFLAGS=
#HPX_LIBRARIES=
#HPX_LDFLAGS=

# NETCDF
NETCDF_ROOT=
NETCDF_CXXFLAGS=
NETCDF_LIBRARIES=
NETCDF_LDFLAGS=

# LAPACK
LAPACK_INTERFACE=
LAPACK_CXXFLAGS=
LAPACK_LIBRARIES=
LAPACK_LDFLAGS=

# DO NOT DELETE
