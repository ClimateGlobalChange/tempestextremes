/////////////////////////////////////////////
///
///          \file blockingUtilities.h
///          \author Marielle Pinheiro
///          \version March 1, 2015
///

#ifndef _BLOCKINGUTILITIES_H_
#define _BLOCKINGUTILITIES_H_

////////////////////////////////////////////////

#include "NetCDFUtilities.h"
#include "netcdfcpp.h"
#include "DataVector.h"
#include "DataMatrix.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"
#include "TimeObj.h"
#include "Announce.h"

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>

//Copied from StitchBlobs
void GetInputFileList(
        const std::string & strInputFileList,
                std::vector<std::string> & vecInputFiles
);

int DayInYear(int nMonth,
              int nDay);

void ParseTimeDouble(
        const std::string & strTimeUnits,
        const std::string & strTimeCalendar,
        double dTime,
        int & nDateYear,
        int & nDateMonth,
        int & nDateDay,
        int & nDateHour
); 



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
  double & lat_res,
  double & lon_res,
  double & dphi,
  double & dlambda,
  double & p_res,
  DataVector<double> & coriolis,
  DataVector<double> & cosphi
);

//Function that calculates PT
void PT_calc(
        NcVar *T, 
        NcVar *pLev, 
        DataMatrix4D<double> &PTMat
);

//Function that calculates relative vorticity
void rVort_calc(
        NcVar *U,
        NcVar *V,
        double dphi,
        double dlambda,
        DataVector<double> cosphi,
        DataMatrix4D<double> & RVMat
);

//Function that calculates PV
void PV_calc(
        NcVar *U,
        NcVar *V,
        DataMatrix4D<double> PTMat,
        DataMatrix4D<double> RVMat,
        NcVar *pVals,
        DataVector<double> coriolis,
        DataVector<double> cosphi,
        double dphi,
        double dlambda,
        double lat_res,
        double lon_res,
        NcVar *PV,
        NcVar *intPV);

void calcDevsPV(bool leap,
              int startAvgIndex,
              NcVar *inIPV,
              NcVar *outDev,
              NcVar *outADev,
              NcVar *outPosIntDev,
              NcVar *avgIPV,
              NcVar *inTime,
              NcVar *avgTime,
              NcVar *lat,
              NcVar *outTime,
              double PVAnom);

void stdDev(DataMatrix3D<double>inDevs,
              int nTime,
              int nLat,
              int nLon,
              DataMatrix<double> & outStdDev);


void calcDevsGH(bool leap,
              int startAvgIndex,
              NcVar *inGH,
              NcVar *outDev,
              NcVar *outADev,
              NcVar *outIntDev,
              NcVar *avgGH,
              NcVar *inTime,
              NcVar *avgTime,
              NcVar *lat,
              NcVar *outTime,
              NcVar *stdDevVar);

//Function that calculates TM blocking index
double GHcheck(double z_0,
          double z_N,
          double z_S,
          double lat_0,
          double lat_N,
          double lat_S,
          std::string hemi );

bool missingValCheck(
  DataMatrix3D<double> fillData,
  int nTime,
  double missingNum
);

bool checkFileLeap(
  std::string StrTimeUnits,
  std::string strCalendar,
  int dateYear,
  int dateMonth,
  int dateDay,
  int dateHour
);


#endif
