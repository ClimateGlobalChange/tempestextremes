//////////////////////////////////
///
///    \file interpolate.h
///    \author Marielle Pinheiro
///    \version March 24, 2015

#ifndef _INTERP_H_
#define _INTERP_H_

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

void interp_util(NcFile & readin,
                 const std::string & strname_2d,
                 NcFile & ifile_out);

#endif
