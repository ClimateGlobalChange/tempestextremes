#!/bin/sh

#ssh -Y cori.nersc.gov            # Cori
salloc -q interactive -N 1 -t 02:00:00 -C haswell
module load arm-forge
module load cray-netcdf
make -j