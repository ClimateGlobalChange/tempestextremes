//////////////////////////////////
///
///    \file interp_z500.h
///    \author Marielle Pinheiro
///    \version November 15, 2015

#ifndef _INTERP_Z_
#define _INTERP_Z_

/////////////////////////////////
#include "BlockingUtilities.h"
#include "NetCDFUtilities.h"
#include "netcdfcpp.h"
#include "DataVector.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"

#include <cstdlib>
#include <cmath>
#include <cstring>

void interp_1lev(NcVar *var,
                     NcVar *hyam,
                     NcVar *hybm,
                     NcVar *ps,
                     double plev,
                     NcVar *NewVar
);

void interp_z500(NcFile & readin,
                 const std::string & strname_2d,
                 NcFile & ifile_out);

#endif
