/////////////////////////////////////////////
///
///          \file CLIVAR_block_utilities.h
///          \author Marielle Pinheiro
///          \version March 1, 2015
///

#ifndef _CLIVARBLOCKUTIL_H_
#define _CLIVARBLOCKUTIL_H_

////////////////////////////////////////////////

#include "NetCDFUtilities.h"
#include "netcdfcpp.h"
#include "DataVector.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"

#include <cstdlib>
#include <cmath>
#include <cstring>

//Function to interpolate variables from hybrid levels to pressure levels
void interpolate_lev(NcVar *var,
                     NcVar *hyam,
                     NcVar *hybm,
                     NcVar *ps,
                     NcVar *pLev,
                     NcVar *NewVar
);

//Function to copy dimension variables to outfile
void copy_dim_var(
        NcVar *inVar,
        NcVar *outVar
);

//Function that calculates dlat, don, etc for PV calculation
void pv_vars_calc(
  NcVar *lat,
  NcVar *lon,
  NcVar *plev,
  double lat_res,
  double lon_res,
  double p_res,
  DataVector<double> coriolis,
  DataVector<double> cosphi
);

//Function that calculates PT
void PT_calc(
        NcVar *T, 
        NcVar *pLev, 
        NcVar *PT
);

//Function that calculates relative vorticity
void rVort_calc(
        NcVar *U,
        NcVar *V,
        double dphi,
        double dlambda,
        DataVector<double> cosphi,
        NcVar *vorticity
);

//Function that calculates PV
void PV_calc(
        NcVar *U,
        NcVar *V,
        NcVar *PT,
        NcVar *rVort,
        DataVector<double> coriolis,
        double dphi,
        double dlambda,
        double dp,
        NcVar *PV);

#endif
