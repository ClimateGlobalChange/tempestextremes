#!/bin/bash
module load cray-hdf5
module load cray-netcdf
cmake -D CMAKE_CXX_COMPILER=CC .
make all
make install
make clean
