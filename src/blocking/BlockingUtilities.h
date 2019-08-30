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
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include "DataArray4D.h"
#include "TimeObj.h"
#include "Announce.h"

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>

//////////////////////////////////////
//    SECTION: FILE OPERATIONS      //
//////////////////////////////////////

//Copied from StitchBlobs: This function
//takes an input text file list and creates
//a character string vector of the file names
void GetInputFileList(
        const std::string & strInputFileList,
                std::vector<std::string> & vecInputFiles
);

				
//////////////////////////////////////
//    SECTION: TIME OPERATIONS      //
//////////////////////////////////////
				
//Used to determine if the input time value
//is a leap day or not. Returns boolean for
//that particular time value
bool checkFileLeap(
  std::string StrTimeUnits,
  std::string strCalendar,
  int dateYear,
  int dateMonth,
  int dateDay,
  int dateHour,
  double timeVal
);

//Function to get number of day in year
//Returns an integer value in the range 1-365
int DayInYear(int nMonth,
              int nDay,
              std::string strCalendar);

//Copied from DetectCyclones: This function
//takes a time value (in units of hours or 
//days since reference date) and returns 4
//integer values: year, month, day, and hour
void ParseTimeDouble(
        const std::string & strTimeUnits,
        const std::string & strTimeCalendar,
        double dTime,
        int & nDateYear,
        int & nDateMonth,
        int & nDateDay,
        int & nDateHour
); 

bool sequentialFiles(
  int prevYear,
  int prevMonth,
  int prevDay,
  int prevHour,
  int nextYear,
  int nextMonth,
  int nextDay,
  int nextHour,
  std::string strCalendar
);

//Takes the last time value of the previous file
//and the first time value of the current file
//and returns the difference in the amount of time 
//between these two values. Used in BlockingAvg to 
//check if this value exceeds the time axis resolution
double tBetweenFiles(
	std::string strTimeUnits,
	double nextStartTime,
	double prevEndTime
);

//////////////////////////////////////
//    SECTION: AXIS OPERATIONS      //
//////////////////////////////////////

//Takes an input file's axis variable and copies
//relevant attributes to output file's axis variable
void copy_dim_var(
        NcVar *inVar,
        NcVar *outVar
);

//Function to interpolate variables from hybrid levels to pressure levels
//Takes input variable and reference variables hyam,hybm and returns
//variable with new pressure level axis (specified by pLev)
void interpolate_lev(NcVar *var,
                     NcVar *hyam,
                     NcVar *hybm,
                     NcVar *ps,
                     NcVar *pLev,
                     NcVar *NewVar
);

//Function that takes an input variable (with pressure axis as vertical)
//and averages variable along the pressure dimension from 150-500 hPa. 
//Returns variable with dimensions [time, lat, lon]
void VarPressureAvg(
    NcVar *invar,
    NcVar * pVals,
    NcVar * outvar
);

/////////////////////////////////////////////////////////
//    SECTION: INTERMEDIATE VARIABLE CALCULATIONS      //
/////////////////////////////////////////////////////////

//Takes temperature and pressure variables and returns
//a 4D data matrix (time, lev, lat, lon) with potential 
//temperature values
void PT_calc(
        int nPlev,
        int nLat,
        int nLon,
        DataArray3D<double> TMat, 
        NcVar *pLev, 
        DataArray3D<double> &PTMat
);

//Replaces a missing value with one interpolated in the longitudinal direction
double replaceMissingFloat(int currA,
                           int currB,
                           int currP,
                           double valThresh,
                           DataArray3D<double> VarMat,
                           int aLen,
                           int bLen
);
double replaceMissingFloat2D(int currA,
                           int currB,
                           double valThresh,
                           DataArray2D<double> VarMat,
                           int aLen,
                           int bLen
);		
//Used in BlockingPV. Input lat, lon, and pressure variables
//and returns the variables necessary to calculate PV (dlat,
//dlon, vector of coriolis parameter values,etc)
void pv_vars_calc(
  NcVar *lat,
  NcVar *lon,
  NcVar *plev,
  double & lat_res,
  double & lon_res,
  double & dphi,
  double & dlambda,
  double & p_res,
  DataArray1D<double> & coriolis,
  DataArray1D<double> & cosphi
);


//Takes wind variables, phi/lambda resolution and coriolis values
//and returns a 4D data matrix (time, lev, lat, lon) with relative 
//vorticity values
void rVort_calc(
        int nPlev,
        int nLat,
        int nLon,
        DataArray3D<double>UMat,
        DataArray3D<double>VMat,
        double dphi,
        double dlambda,
        DataArray1D<double> cosphi,
        DataArray3D<double> & RVMat
);

//////////////////////////////////////
//    SECTION: VARIABLE CHECKS      //
//////////////////////////////////////


double GHcheck(double z_0,
          double z_N,
          double z_S,
          double lat_0,
          double lat_N,
          double lat_S,
          std::string hemi );




bool missingValCheck(
  DataArray3D<double> fillData,
  int nTime,
  double missingNum
);


void MissingFill(
  double missingValue,
  double tRes,
  double contCheck,
  int nLat,
  int nLon,
  int ArrLen,
  int & currArrIndex,
  int & dateIndex,
  DataArray3D<double> & currFillData
);

//////////////////////////////////////////////////
//    SECTION: FINAL VARIABLE CALCULATIONS      //
//////////////////////////////////////////////////

//Takes  variables for wind, potential temperature,
//relative vorticity, etc and outputs both 3D 
//(lev, lat, lon) and 2D (lat, lon) vertically 
//averaged potential vorticity variables per time slice

void PV_calc(
        int nPlev,
        int nLat,
        int nLon,
        DataArray3D<double>UMat,
        DataArray3D<double>VMat,
        DataArray3D<double> PTMat,
        DataArray3D<double> RVMat,
        DataArray1D<double>pVec,
        DataArray1D<double> coriolis,
        DataArray1D<double>cosphi,
        double dphi,
        double dlambda,
        double lat_res,
        double lon_res,
        DataArray3D<double> & PVMat
);


void IPV_calc(
       int nPlev,
       int nLat,
       int nLon,
       double lat_res,
       DataArray1D<double> pVec,
       DataArray3D<double> PVMat,
       DataArray2D<double> & IPVMat
);

////////////////////////////////////////////////////
//    SECTION: VARIABLE ANOMALY CALCULATIONS      //
////////////////////////////////////////////////////

//Takes instantaneous variable (PV or Z) and long term 
//daily average and outputs 3 variables: instantaneous 
//anomalies, anomalies with 2-day smoothing, and a normalized
//anomaly (all values below threshold or wrong sign set to 0)
/*void calcDevs(bool isPV,
              std::string ZtoGH,
              std::string is4D,
              int pIndex,
              int nSteps,
              int nOutTime,
              std::string strTimeUnits,
              std::string strCalendar,
              NcVar *inIPV,
              NcVar *outDev,
              NcVar *outADev,
              NcVar *avgIPV,
              NcVar *inTime,
              NcVar *avgTime,
              NcVar *lat);
*/

void calcDevs(bool latNorm,
              std::string ZtoGH,
              std::string is4D,
              int pIndex,
              int nSteps,
              int nOutTime,
              std::string strTimeUnits,
              std::string strCalendar,
              NcVar *inIPV,
              NcVar *outDev,
              NcVar *avgIPV,
              NcVar *inTime,
              NcVar *avgTime,
              NcVar *lat,
              double missingNo);
							
/*void calcSmoothedDevs( NcVar *outADev,
											int nTime,
											int nLat,
											int nLon,
											int nSteps,
						          int & currMatIndex,
						          int & currTindex,
						          DataArray3D<double> &twoDayMat);
*/
void calcNormalizedDevs(bool isPV,
                       NcVar * inDev,
                       NcVar * outPosIntDev,
                       NcVar * lat,
                       NcVar * inTime,
                       std::string strTimeUnits,
                       std::string strCalendar,
                       DataArray3D<double>threshMat,
                       double minThresh);
/*void stdDev(DataArray3D<double>inDevs,
              int nTime,
              int nLat,
              int nLon,
              DataArray2D<double> & outStdDev);


void calcDevsGH(bool leap,
              double GHAnom,
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
*/
//Function that calculates TM blocking index



void PV_calc2(
        int nPlev,
        int nLat,
        int nLon,
        DataArray3D<double>UMat,
        DataArray3D<double>VMat,
        DataArray3D<double> PTMat,
        DataArray3D<double> RVMat,
        DataArray1D<double>pVec,
        DataArray1D<double> coriolis,
        DataArray1D<double>cosphi,
        double dphi,
        double dlambda,
        double lat_res,
        double lon_res,
        DataArray3D<double> &PVMat,
  DataArray3D<double> &dpt_dp,
  DataArray3D<double> &du_dp,
  DataArray3D<double> &dv_dp,
  DataArray3D<double> &dpt_dphi,
  DataArray3D<double> &dpt_dl
);

#endif
