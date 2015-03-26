//////////////////////////////////
///
///    \file interpolate.h
///    \author Marielle Pinheiro
///    \version March 24, 2015

#ifndef _INTERP_H_
#define _INTERP_H_

/////////////////////////////////
#include "CLIVAR_block_utilities.h"
#include "NetCDFUtilities.h"
#include "netcdfcpp.h"
#include "DataVector.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"

#include <cstdlib>
#include <cmath>
#include <cstring>

void interp_util(NcFile readin,
                 std::string strname_2d,
                 std::string interp_out);

#endif
