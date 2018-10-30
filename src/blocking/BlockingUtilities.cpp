/////////////////////////////////////////////
///
///          \file blockingUtilities.cpp
///
///          \author Marielle Pinheiro
///          
///          \version June 1, 2015
///

#include "BlockingUtilities.h"
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

//////////////////////////////////////
//    SECTION: FILE OPERATIONS      //
//////////////////////////////////////

//Copied from StitchBlobs: This function
//takes an input text file list and creates
//a character string vector of the file names+
void GetInputFileList(
        const std::string & strInputFileList,
                std::vector<std::string> & vecInputFiles
        ) {
                FILE * fp = fopen(strInputFileList.c_str(), "r");

                char szBuffer[1024];
                for (;;) {
                        fgets(szBuffer, 1024, fp);

                        if (feof(fp)) {
                                break;
                        }

                        // Remove end-of-line characters
                        for (;;) {
                                int nLen = strlen(szBuffer);
                                if ((szBuffer[nLen-1] == '\n') ||
                                        (szBuffer[nLen-1] == '\r') ||
                                        (szBuffer[nLen-1] == ' ')
                                ) {
                                        szBuffer[nLen-1] = '\0';
                                        continue;
                                }
                                break;
                        }

                        vecInputFiles.push_back(szBuffer);
        }

        if (vecInputFiles.size() == 0) {
                _EXCEPTION1("No files found in file \"%s\"", strInputFileList.c_str());
        }

        fclose(fp);
}

//////////////////////////////////////
//    SECTION: TIME OPERATIONS      //
//////////////////////////////////////

//Used to determine if the input time value
//is a leap day or not. Returns boolean for
//that particular time value
bool checkFileLeap(
  std::string strTimeUnits,
  std::string strCalendar,
  int dateYear,
  int dateMonth,
  int dateDay,
  int dateHour,
  double timeVal
){

  bool leap = false;

  int leapYear=0;
  int leapMonth=0;
  int leapDay=0;
  int leapHour=0;

  if (strCalendar!="noleap" && dateMonth<=2){
    //Check whether file contains a Feb 29

    ParseTimeDouble(strTimeUnits, strCalendar, timeVal, leapYear,\
      leapMonth, leapDay, leapHour);

    if ((leapMonth==2 && leapDay==29) || (dateMonth==2&&leapMonth==3)){
      //Check when parsing the indices
      leap = true;
    }
  }
  return leap;
}

//Function to get number of day in year
//Returns an integer value in the range 1-365
int DayInYear(int nMonth, int nDay){
  int day=0;
  if (nMonth>1){
    for (int x=1; x<nMonth; x++){
      if (x==2){
        day += 28;
      }
      else if (x==4 || x==6 || x==9 || x==11){
        day +=30;
      }
      else{
        day +=31;
      }
    }
  }
  day+=nDay;
  return(day);
}

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
) {
	// Get calendar type
	Time::CalendarType cal;
	if ((strTimeCalendar.length() >= 6) &&
		(strncmp(strTimeCalendar.c_str(), "noleap", 6) == 0)
	) {
		cal = Time::CalendarNoLeap;

	} else if (
		(strTimeCalendar.length() >= 8) &&
		(strncmp(strTimeCalendar.c_str(), "standard", 8) == 0)
	) {
		cal = Time::CalendarStandard;

	} else if (
                (strTimeCalendar.length() >=9) &&
                (strncmp(strTimeCalendar.c_str(), "gregorian",9)==0)
        ) {
                cal = Time::CalendarStandard;
        } else {
		_EXCEPTION1("Unknown calendar type \"%s\"", strTimeCalendar.c_str());
	}
	// Time format is "days since ..."
	if ((strTimeUnits.length() >= 11) &&
	    (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)
	) {
		std::string strSubStr = strTimeUnits.substr(11);
		Time time(cal);
		time.FromFormattedString(strSubStr);

		int nDays = static_cast<int>(dTime);
		time.AddDays(nDays);

		int nSeconds = static_cast<int>(fmod(dTime, 1.0) * 86400.0);
		time.AddSeconds(nSeconds);

		nDateYear = time.GetYear();
		nDateMonth = time.GetMonth();
		nDateDay = time.GetDay();
		nDateHour = time.GetSecond() / 3600;

	// Time format is "hours since ..."
	} else if (
	    (strTimeUnits.length() >= 12) &&
	    (strncmp(strTimeUnits.c_str(), "hours since ", 12) == 0)
	) {
		std::string strSubStr = strTimeUnits.substr(12);
		Time time(cal);
		time.FromFormattedString(strSubStr);
		time.AddHours(static_cast<int>(dTime));

		nDateYear = time.GetYear();
		nDateMonth = time.GetMonth();
		nDateDay = time.GetDay();
		nDateHour = time.GetSecond() / 3600;

	// Time format is "minutes since ..."
	} else if (
	    (strTimeUnits.length() >= 14) &&
	    (strncmp(strTimeUnits.c_str(), "minutes since ", 14) == 0)
	) {
		std::string strSubStr = strTimeUnits.substr(14);
		Time time(cal);
		time.FromFormattedString(strSubStr);
		time.AddMinutes(static_cast<int>(dTime));

		nDateYear = time.GetYear();
		nDateMonth = time.GetMonth();
		nDateDay = time.GetDay();
		nDateHour = time.GetSecond() / 3600;

	} else {
		_EXCEPTIONT("Unknown \"time::units\" format");
	}
}

bool sequentialFiles(
  int prevYear,
  int prevMonth,
  int prevDay,
  int prevHour,
  int nextYear,
  int nextMonth,
  int nextDay,
  int nextHour
){
  bool isSequential = false;
  //Day in year of each date?
  int dayPrev = DayInYear(prevMonth,prevDay);
  int dayNext = DayInYear(nextMonth,nextDay);
  //Difference between years?
  int yearDiff = nextYear-prevYear;
  //If yearDiff is 1, is it going from December to January?
  int dayDiff = dayNext-dayPrev;
  if (yearDiff <=1.){
    //Is the difference 365-1?
    if (std::fabs(dayDiff)>=364){
      isSequential = true;
    }
    if (std::fabs(dayDiff)<=1){
      isSequential = true;
    }
  }
  return isSequential;
}

//Takes the last time value of the previous file
//and the first time value of the current file
//and returns the difference in the amount of time 
//between these two values. Used in BlockingAvg to 
//check if this value exceeds the time axis resolution
double tBetweenFiles(
  std::string strTimeUnits,
  double nextStartTime,
  double prevEndTime
){
  double contCheck;
  if ((strTimeUnits.length() >= 11) && \
    (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)){
    contCheck = std::fabs(nextStartTime-prevEndTime);
  }
  else{
    contCheck = std::fabs(nextStartTime-prevEndTime)/24.0;
  }
  return(contCheck);
}

//////////////////////////////////////
//    SECTION: AXIS OPERATIONS      //
//////////////////////////////////////

//Takes an input file's axis variable and copies
//relevant attributes to output file's axis variable
void copy_dim_var(
	NcVar *inVar, 
	NcVar *outVar
	){
  //Get necessary information from infile
  int varLen = inVar->get_dim(0)->size();
  DataVector<double> inVec(varLen);
  inVar->set_cur((long) 0);
  inVar->get(&(inVec[0]), varLen);
  //Copy data to new outgoing variable
  outVar->set_cur((long) 0);
  outVar->put(&(inVec[0]), varLen);
  //Copy other attributes
  CopyNcVarAttributes(inVar, outVar);
}

//Function to interpolate variables from hybrid levels to pressure levels
//Takes input variable and reference variables hyam,hybm and returns
//variable with new pressure level axis (specified by pLev)
void interpolate_lev(NcVar *var, 
                     NcVar *hyam, 
                     NcVar *hybm, 
                     NcVar *ps, 
                     NcVar *pLev,
                     NcVar *NewVar 
){

  int nTime,nLev,nLat,nLon,npLev;
  double A1,A2,B1,B2,p1,p2,weight;

  //Array dimensions
  nTime = var->get_dim(0)->size();
  nLev = var->get_dim(1)->size();
  nLat = var->get_dim(2)->size();
  nLon = var->get_dim(3)->size();
  npLev = pLev->get_dim(0)->size();


  //hybrid coefficient A
  DataVector<double> vecHyam(nLev);
  hyam->set_cur((long) 0);
  hyam->get(&(vecHyam[0]), nLev);

  //hybrid coefficient B
  DataVector<double> vecHybm(nLev);
  hybm->set_cur((long) 0);
  hybm->get(&(vecHybm[0]), nLev);

  //Pressure levels
  DataVector<double> vecpLev(npLev);
  pLev->set_cur((long) 0);
  pLev->get(&(vecpLev[0]), npLev);
 
  //Matrix to store PS
  DataMatrix <double> matPS(nLat, nLon);

  //Matrix to store input variable data
  DataMatrix3D<double> matVar(nLev, nLat, nLon);
  
  //Matrix to store output variable data
  DataMatrix3D<double> matVarOut(npLev, nLat, nLon);
  //std::cout<<"within interpolate_lev: about to interpolate"<<std::endl; 
  //Loop over input data and interpolate to output var
  for (int t=0; t<nTime; t++){
    ps->set_cur(t, 0, 0);
    ps->get(&(matPS[0][0]), 1, nLat, nLon);
    var->set_cur(t, 0, 0, 0);
    var->get(&(matVar[0][0][0]), 1, nLev, nLat, nLon);
    for (int l=0; l<(nLev-1); l++){
      A1 = vecHyam[l];
      B1 = vecHybm[l];
      A2 = vecHyam[l+1];
      B2 = vecHybm[l+1];
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
	  int interp_check = 0;
          p1 = 100000.0 * A1 + matPS[a][b] * B1;
          p2 = 100000.0 * A2 + matPS[a][b] * B2;
          for (int p=0; p<npLev; p++){
            if (p1<vecpLev[p] && p2>=vecpLev[p]){
	      interp_check = 1;
              weight = ((vecpLev[p]-p1)/(p2-p1));
              matVarOut[p][a][b] = weight*matVar[l+1][a][b]
                               + (1.0-weight)*matVar[l][a][b];
            }
          }
        }
      }
    }
    NewVar->set_cur(t, 0, 0, 0);
    NewVar->put(&(matVarOut[0][0][0]),1, npLev, nLat, nLon);

  }
  CopyNcVarAttributes(var, NewVar);
}

//Function that takes an input variable (with pressure axis as vertical)
//and averages variable along the pressure dimension from 150-500 hPa. 
//Returns variable with dimensions [time, lat, lon]
void VarPressureAvg(
	NcVar * invar,
	NcVar * pVals,
	NcVar * outvar
){
	int nTime,nLat,nLon,nPlev;
	nTime = invar->get_dim(0)->size();
	nPlev = invar->get_dim(1)->size();
	nLat = invar->get_dim(2)->size();
	nLon = invar->get_dim(3)->size();

	//Input into Matrix
	DataMatrix4D<double> inMat(nTime,nPlev,nLat,nLon);
	invar->set_cur(0,0,0,0);
	invar->get(&(inMat[0][0][0][0]),nTime,nPlev,nLat,nLon);
	
    //Pressure axis values
    DataVector<double> pVec(nPlev);
    pVals->set_cur((long) 0);
    pVals->get(&(pVec[0]), nPlev);

	
	//Integrate values over upper troposphere
    int pos_top;
    int pos_bot;

    for (int x=0; x<nPlev; x++){
      if (std::fabs(pVec[x]-15000.0)<0.0001){
        pos_top = x;
      }
      if (std::fabs(pVec[x]-50000.0)<0.0001){
        pos_bot = x;
      }
    }
   //if level axis starts at ground, reverse positions
    if (pos_top>pos_bot){
      int temp = pos_bot;
      pos_bot = pos_top;
      pos_top = temp;
    }

    DataMatrix3D<double> outMat(nTime, nLat, nLon);  
    double bot,mid,top;
    double modLevLen = pos_bot-pos_top;
    double invLevLen = 1.0/(2.0*modLevLen);

    //Calculate integration parts
    for (int t=0; t<nTime; t++){
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          top = inMat[t][pos_top][a][b];
          bot = inMat[t][pos_bot][a][b];
          mid = 0.0;
          for (int p=(pos_top+1); p<pos_bot; p++){
            mid+=2.0*inMat[t][p][a][b];
          }
          outMat[t][a][b] = (top+bot+mid)*invLevLen;
        }
      }
    }
 
    outvar->set_cur(0,0,0);
    outvar->put(&(outMat[0][0][0]),nTime,nLat,nLon);
    std::cout<<"Finished integrating variable."<<std::endl;
}

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
	DataMatrix3D<double>TMat, 
	NcVar *pLev, 
	DataMatrix3D<double> &PTMat
){

  double pFrac;
  double exp = 287.0/1004.5;

  //INPUTS
  DataVector<double> pVec(nPlev);

  pLev->set_cur((long) 0);
  pLev->get(&(pVec[0]), nPlev);
  

    int currA = 0;
    int prevA = 0;
    int nextA = 0;
    int currB = 0;
    int prevB = 0;
    int nextB = 0;
    double currT;
  //OUTPUT: PT
//  DataMatrix4D<double> PTMat(nTime, nPlev, nLat, nLon);
    for (int p=0; p<nPlev; p++){
      pFrac = 100000.0/pVec[p];
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          currT = TMat[p][a][b];
          if (std::fabs(currT)>1000.0){
            currT = replaceMissingFloat(a,b,p,1000.0,TMat,nLat,nLon);
          }
          PTMat[p][a][b] = currT*std::pow(pFrac, exp);
        }
      }
    }
}


double replaceMissingFloat(int currA,
                           int currB,
                           int currP,
                           double valThresh,
                           DataMatrix3D<double> VarMat,
                           int aLen,
                           int bLen
){
  double currVal = VarMat[currP][currA][currB];
  //Point left
  int leftB = currB-1;
  double leftVal = VarMat[currP][currA][leftB];
  while (std::fabs(leftVal) >= std::fabs(valThresh)){
    leftB-=1;
    if (leftB < 1){
      break;
    }
    leftVal = VarMat[currP][currA][leftB];
  }
  //Point right
  int rightB = currB+1;
  double rightVal = VarMat[currP][currA][rightB];
  while (std::fabs(rightVal)>= std::fabs(valThresh)){
    rightB+=1;
    if (rightB > (bLen-1)){
      break;
    }
    rightVal = VarMat[currP][currA][rightB];
  }

  //Interpolate!
  double newCurrVal = ((double)(rightB-currB)/(double)(rightB-leftB))*leftVal + \
     ((double)(currB-leftB)/(double)(rightB-leftB))*rightVal;
  return(newCurrVal);
}

double replaceMissingFloat2D(int currA,
                           int currB,
                           double valThresh,
                           DataMatrix<double> VarMat,
                           int aLen,
                           int bLen
){
  double currVal = VarMat[currA][currB];
  //Point left
  int leftB = currB-1;
  double leftVal = VarMat[currA][leftB];
  while (std::fabs(leftVal) >=std::fabs(valThresh)){
    leftB-=1;
    if (leftB < 1){
      break;
    }
    leftVal = VarMat[currA][leftB];
  }
  //Point right
  int rightB = currB+1;
  double rightVal = VarMat[currA][rightB];
  while (std::fabs(rightVal)>=std::fabs(valThresh)){
    rightB+=1;
    if (rightB > (bLen-1)){
      break;
    }
    rightVal = VarMat[currA][rightB];
  }

  //Interpolate!
  double newCurrVal = ((double)(rightB-currB)/(double)(rightB-leftB))*leftVal + \
     ((double)(currB-leftB)/(double)(rightB-leftB))*rightVal;
  return(newCurrVal);
}

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
  DataVector<double> & coriolis,
  DataVector<double> & cosphi
){ 
  double radius = 6371000.0;
  double pi = 4.0*std::atan(1.0);
  double inv_radian = pi/180.0;
  double sigma = 7.2921*std::pow(10.0, -5.0);
  int nLat,nLon,np;
  double sinphi;

  nLat = lat->get_dim(0)->size();
  DataVector<double> latVec(nLat);
  lat->set_cur((long) 0);
  lat->get(&(latVec[0]), nLat);

  nLon = lon->get_dim(0)->size();
  DataVector<double> lonVec(nLon);
  lon->set_cur((long) 0);
  lon->get(&(lonVec[0]), nLon);

  np = plev->get_dim(0)->size();
  DataVector<double> pVec(np);
  plev->set_cur((long) 0);
  plev->get(&(pVec[0]), np);
  
 
  lat_res = latVec[2]-latVec[1];
  lon_res = lonVec[2]-lonVec[1];
  dphi = lat_res*inv_radian;
  dlambda = lon_res*inv_radian;
  p_res = pVec[2]-pVec[1];
  for (int i=0; i<nLat; i++){
    sinphi = std::sin(latVec[i]*inv_radian);
    coriolis[i] = 2.0 * sigma * sinphi;
    cosphi[i] = std::cos(latVec[i]*inv_radian);
  }
}

//Takes wind variables, phi/lambda resolution and coriolis values
//and returns a 4D data matrix (time, lev, lat, lon) with relative 
//vorticity values
void rVort_calc(
        int nPlev,
        int nLat,
        int nLon,
	DataMatrix3D<double>UMat,
	DataMatrix3D<double>VMat,
	double dphi,
	double dlambda,
	DataVector<double> cosphi,
	DataMatrix3D<double> & RVMat
){

  double radius = 6371000.0;
  double invDlambda = 1.0/(2.0*dlambda);
  double invDphi = 1.0/(2.0*dphi);
  double coef;

//OUTPUT: Partial derivatives
//U WRT PHI
  double U1,U2,U3,V1,V2,V3;

  DataMatrix3D<double> dUdphi(nPlev, nLat, nLon);
    for (int p=0; p<nPlev; p++){
      for (int b=0; b<nLon; b++){
        U1=UMat[p][2][b];
        if (std::fabs(U1)>1000.0){
          U1=replaceMissingFloat(2,b,p,1000.0,UMat,nLat,nLon);
        }
        U2=UMat[p][1][b];
        if (std::fabs(U2)>1000.0){
          U2=replaceMissingFloat(1,b,p,1000.0,UMat,nLat,nLon);
        }
        U3=UMat[p][0][b];
        if (std::fabs(U3)>1000.0){
          U3=replaceMissingFloat(0,b,p,1000.0,UMat,nLat,nLon);
        }

        dUdphi[p][0][b]=(-U1*cosphi[2]+4.0*U2*cosphi[1]\
          -3.0*U3*cosphi[0])*invDphi;


        U1=UMat[p][nLat-1][b];
        if (std::fabs(U1)>1000.0){
          U1=replaceMissingFloat((nLat-1),b,p,1000.0,UMat,nLat,nLon);
        }
        U2=UMat[p][nLat-2][b];
        if (std::fabs(U2)>1000.0){
          U2=replaceMissingFloat((nLat-2),b,p,1000.0,UMat,nLat,nLon);
        }
        U3=UMat[p][nLat-3][b];
        if (std::fabs(U3)>1000.0){
          U3=replaceMissingFloat((nLat-3),b,p,1000.0,UMat,nLat,nLon);
        }
        dUdphi[p][nLat-1][b]=(3.0*U1*cosphi[nLat-1]-4.0*U2*cosphi[nLat-2]\
          +U3*cosphi[nLat-3])*invDphi;
        for (int a=1; a<(nLat-1); a++){
          U1=UMat[p][a+1][b];
          if (std::fabs(U1)>1000.0){
            U1=replaceMissingFloat((a+1),b,p,1000.0,UMat,nLat,nLon);
          }
          U2=UMat[p][a-1][b];
          if (std::fabs(U2)>1000.0){
            U2=replaceMissingFloat((a-1),b,p,1000.0,UMat,nLat,nLon);
          }
          dUdphi[p][a][b]=(U1*cosphi[a+1]\
           -U2*cosphi[a-1])*invDphi;
        }
      }
    }

  //V WRT LAMBDA
  DataMatrix3D<double> dVdl(nPlev, nLat, nLon);
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        V1=VMat[p][a][1];
        if (std::fabs(V1)>1000.0){
          V1=replaceMissingFloat(a,1,p,1000.0,VMat,nLat,nLon);
        }
        V2=VMat[p][a][nLon-1];
        if (std::fabs(V2)>1000.0){
          V2=replaceMissingFloat(a,(nLon-1),p,1000.0,VMat,nLat,nLon);
        }
        dVdl[p][a][0]=(V1-V2)*invDlambda;
        V1=VMat[p][a][0];
        if (std::fabs(V1)>1000.0){
          V1=replaceMissingFloat(a,0,p,1000.0,VMat,nLat,nLon);
        }
        V2=VMat[p][a][nLon-2];
        if (std::fabs(V2)>1000.0){
          V2=replaceMissingFloat(a,(nLon-2),p,1000.0,VMat,nLat,nLon);
        }
        dVdl[p][a][nLon-1]=(V1-V2)*invDlambda;
        for (int b=1; b<(nLon-1); b++){
          V1=VMat[p][a][b+1];
          if (std::fabs(V1)>1000.0){
            V1=replaceMissingFloat(a,(b+1),p,1000.0,VMat,nLat,nLon);
          }
          V2=VMat[p][a][b-1];
          if (std::fabs(V2)>1000.0){
            V2=replaceMissingFloat(a,(b-1),p,1000.0,VMat,nLat,nLon);
          }

          dVdl[p][a][b]=(V1-V2)*invDlambda;
        }
      }
    }

    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        coef = 1.0/(radius*cosphi[a]);
        for (int b=0; b<nLon; b++){
          RVMat[p][a][b] = coef*(dVdl[p][a][b]-dUdphi[p][a][b]);
        }
      }
    }

  //Lat end cases: set to 0 because of pole singularities causing error
    for (int p=0; p<nPlev; p++){
      for (int b=0; b<nLon; b++){
        RVMat[p][0][b] = 0.0;
        RVMat[p][nLat-1][b] = 0.0;
      }
    }

//  std::cout<<"Finished calculating relative vorticity."<<std::endl;
}

//////////////////////////////////////
//    SECTION: VARIABLE CHECKS      //
//////////////////////////////////////


double GHcheck(double z_0,
          double z_N,
          double z_S,
          double lat_0,
          double lat_N,
          double lat_S,
          std::string hemi ){
  double GHGS;
  double GHGN;
  double gval;
  //NH: GHGS>0, GHGN<-10m/deg lat
  if (hemi =="N"){
    GHGS = (z_0-z_S)/(lat_0-lat_S);
    GHGN = (z_N-z_0)/(lat_N-lat_0);

    if ((GHGS>0) && (GHGN < -10)){
      gval=1.;
    }
    else gval=0.;
  //SH: GHGN>0,GHGS<-10m/deg lat
  }else if (hemi=="S"){
    GHGS =-(z_S-z_0)/(lat_S-lat_0);
    GHGN =-(z_0-z_N)/(lat_0-lat_N);

    if ((GHGN>0) && (GHGS<-10)){
      gval=1.;
    }
    else gval=0.;
  }
  else std::cout << "Error: invalid hemisphere specified."<<std::endl;
  return(gval);
}




bool missingValCheck(
  DataMatrix3D<double> fillData,
  int nTime,
  double missingNum
){
    bool isMissing = false;
    for (int t=0; t<nTime; t++){
      if (fillData[t][2][2] == missingNum){
        isMissing = true;
        break;
      }
    }
    if (isMissing == false) std::cout << "Did not find missing values."<<std::endl;
    return isMissing;
}

void MissingFill(
  double missingValue,
  double tRes,
  double contCheck,
  int nLat,
  int nLon,
  int arrLen,
  int & currArrIndex,
  int & dateIndex,
  DataMatrix3D<double> & currFillData
){

  int nFill = contCheck/tRes-1;
  int nDaysSkip = int(contCheck-tRes);
  //std::cout<<"ContCheck is " << contCheck << " and tRes is "<<\
   tRes <<" and nFill is "<<nFill << std::endl;
  for (int n=0; n<nFill; n++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        currFillData[currArrIndex][a][b] = missingValue;
      } 
    }
    currArrIndex +=1;
    if (currArrIndex >= arrLen){
      currArrIndex -= arrLen;
    }
  }
  dateIndex += nDaysSkip;
  if (dateIndex >= 365){
    dateIndex-=365;
  }
}


//////////////////////////////////////////////////
//    SECTION: FINAL VARIABLE CALCULATIONS      //
//////////////////////////////////////////////////

//Takes 4D variables for wind, potential temperature,
//relative vorticity, etc and outputs both 4D (time,
//lev, lat, lon) and 3D (time, lat, lon) vertically 
//averaged potential vorticity variables
void PV_calc(
        int nPlev,
        int nLat,
        int nLon,
	DataMatrix3D<double>UMat,
	DataMatrix3D<double>VMat,
	DataMatrix3D<double> PTMat,
	DataMatrix3D<double> RVMat,
        DataVector<double>pVec,	
	DataVector<double> coriolis,
        DataVector<double>cosphi,
	double dphi,
	double dlambda,
        double lat_res,
        double lon_res,
	DataMatrix3D<double> &PVMat
){
  double invdp,invdp1,invdp2;
  double invdphi= 1.0/(2.0*dphi);
  double invdlambda = 1.0/(2.0*dlambda);
  double radius = 6371000.0;
  double coef1,coef2,corvar;
  double U1,U2,U3,V1,V2,V3;


  //Matrices for the partials
  //PT, U, V WRT P
  DataMatrix3D<double> dpt_dp(nPlev, nLat, nLon);
  DataMatrix3D<double> du_dp(nPlev, nLat, nLon);
  DataMatrix3D<double> dv_dp(nPlev, nLat, nLon);

  invdp1 = 1.0/(2.0*(pVec[1]-pVec[0]));
  invdp2 = 1.0/(2.0*(pVec[nPlev-1]-pVec[nPlev-2]));

    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        U1=UMat[2][a][b];
        if (std::fabs(U1)>1000.0){
          U1 = replaceMissingFloat(a,b,2,1000.0,UMat,nLat,nLon);
        }
        U2=UMat[1][a][b];
        if (std::fabs(U2)>1000.0){
          U2 = replaceMissingFloat(a,b,1,1000.0,UMat,nLat,nLon);
        }
        U3=UMat[0][a][b];
        if (std::fabs(U3)>1000.0){
          U3 = replaceMissingFloat(a,b,0,1000.0,UMat,nLat,nLon);
        }
       
        V1=VMat[2][a][b];
        if (std::fabs(V1)>1000.0){
          V1 = replaceMissingFloat(a,b,2,1000.0,VMat,nLat,nLon);
        }
        V2=VMat[1][a][b];
        if (std::fabs(V2)>1000.0){
          V2 = replaceMissingFloat(a,b,1,1000.0,VMat,nLat,nLon);
        }
        V3=VMat[0][a][b];
        if (std::fabs(V3)>1000.0){
          V3 = replaceMissingFloat(a,b,0,1000.0,VMat,nLat,nLon);
        }
      //0 case
        dpt_dp[0][a][b] = (-PTMat[2][a][b]+4.0*PTMat[1][a][b]-3.0*PTMat[0][a][b])*invdp1;
        du_dp[0][a][b] = (-U1+4.0*U2-3.0*U3)*invdp1;
        dv_dp[0][a][b] = (-V1+4.0*V2-3.0*V3)*invdp1;          
      //end case
        U1=UMat[nPlev-1][a][b];
        if (std::fabs(U1)>1000.0){
          U1 = replaceMissingFloat(a,b,(nPlev-1),1000.0,UMat,nLat,nLon);
        }
        U2=UMat[nPlev-2][a][b];
        if (std::fabs(U2)>1000.0){
          U2 = replaceMissingFloat(a,b,(nPlev-2),1000.0,UMat,nLat,nLon);
        }
        U3=UMat[nPlev-3][a][b];
        if (std::fabs(U3)>1000.0){
          U3 = replaceMissingFloat(a,b,(nPlev-3),1000.0,UMat,nLat,nLon);
        }
       
        V1=VMat[nPlev-1][a][b];
        if (std::fabs(V1)>1000.0){
          V1 = replaceMissingFloat(a,b,(nPlev-1),1000.0,VMat,nLat,nLon);
        }
        V2=VMat[nPlev-2][a][b];
        if (std::fabs(V2)>1000.0){
          V2 = replaceMissingFloat(a,b,nPlev-2,1000.0,VMat,nLat,nLon);
        }
        V3=VMat[nPlev-3][a][b];
        if (std::fabs(V3)>1000.0){
          V3 = replaceMissingFloat(a,b,nPlev-3,1000.0,VMat,nLat,nLon);
        }
        dpt_dp[nPlev-1][a][b] = (3.0*PTMat[nPlev-1][a][b]-4.0*PTMat[nPlev-2][a][b]\
          +PTMat[nPlev-3][a][b])*invdp2;
        du_dp[nPlev-1][a][b] = (3.0*U1-4.0*U2+U3)*invdp2;
        dv_dp[nPlev-1][a][b] = (3.0*V2-4.0*V2+V3)*invdp2;
        for (int p=1; p<(nPlev-1); p++){
          U1=UMat[p+1][a][b];
          if (std::fabs(U1)>1000.0){
            U1 = replaceMissingFloat(a,b,(p+1),1000.0,UMat,nLat,nLon);
          }
          U2=UMat[p-1][a][b];
          if (std::fabs(U2)>1000.0){
            U2 = replaceMissingFloat(a,b,(p-1),1000.0,UMat,nLat,nLon);
          }
       
          V1=VMat[p+1][a][b];
          if (std::fabs(V1)>1000.0){
            V1 = replaceMissingFloat(a,b,(p-1),1000.0,VMat,nLat,nLon);
          }
          V2=VMat[p-1][a][b];
          if (std::fabs(V2)>1000.0){
            V2 = replaceMissingFloat(a,b,(p-1),1000.0,VMat,nLat,nLon);
          }
          invdp = 1.0/(2.0*(pVec[p+1]-pVec[p]));
          dpt_dp[p][a][b] = (PTMat[p+1][a][b]-PTMat[p-1][a][b])*invdp;
          du_dp[p][a][b] = (U1-U2)*invdp;
          dv_dp[p][a][b] = (V1-V2)*invdp;
        }
      }
    }
  
  //PT WRT PHI
  DataMatrix3D<double> dpt_dphi(nPlev, nLat, nLon);
  //end cases
    for (int p=0; p<nPlev; p++){
      for (int b=0; b<nLon; b++){
        dpt_dphi[p][0][b]=(-PTMat[p][2][b]+4.0*PTMat[p][1][b]\
          -3.0*PTMat[p][0][b])*invdphi;
        dpt_dphi[p][nLat-1][b]=(3.0*PTMat[p][nLat-1][b]-4.0*PTMat[p][nLat-2][b]\
          +PTMat[p][nLat-2][b])*invdphi;
        for (int a=1; a<(nLat-1); a++){
          dpt_dphi[p][a][b]=(PTMat[p][a+1][b]-PTMat[p][a-1][b])*invdphi;
        }
      }
    }

  //PT WRT LAMBDA
  DataMatrix3D<double> dpt_dl(nPlev, nLat, nLon);
  //end cases
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        dpt_dl[p][a][0]=(PTMat[p][a][1]-PTMat[p][a][nLon-1])*invdlambda;
        dpt_dl[p][a][nLon-1]=(PTMat[p][a][nLon-2]-PTMat[p][a][0])*invdlambda;
        for (int b=1; b<(nLon-1); b++){
          dpt_dl[p][a][b]=(PTMat[p][a][b+1]-PTMat[p][a][b-1])*invdlambda;
        }
      }
    }
  coef2 = 1.0/radius;
  //PV Calculation!
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        coef1 = 1.0/(radius*cosphi[a]);
        corvar=coriolis[a];
        for (int b=0; b<nLon; b++){
          PVMat[p][a][b] = 9.80616*(coef1*dv_dp[p][a][b]*dpt_dl[p][a][b]\
           -coef2*du_dp[p][a][b]*dpt_dphi[p][a][b]\
           -(corvar+RVMat[p][a][b])*dpt_dp[p][a][b]);
           if (std::fabs(PVMat[p][a][b])>10){
             std::cout<<"Questionable value at p="<<p<<", a="<<a<<", b="<<b<<std::endl;
           }
        }
      }
    }
} 

void IPV_calc(
       int nPlev,
       int nLat,
       int nLon,
       double lat_res,
       DataVector<double> pVec,
       DataMatrix3D<double> PVMat,
       DataMatrix<double> & IPVMat       
){
  //Integrate PV over upper troposphere
  int pos_top;
  int pos_bot;

  for (int x=0; x<nPlev; x++){
    if (std::fabs(pVec[x]-15000.0)<0.0001){
      pos_top = x;
    }
    if (std::fabs(pVec[x]-50000.0)<0.0001){
      pos_bot = x;
    }
  }
 //if level axis starts at ground, reverse positions
  if (pos_top>pos_bot){
    int temp = pos_bot;
    pos_bot = pos_top;
    pos_top = temp;
  }

  double PVbot;
  double PVmid;
  double PVtop;

//Eliminate polar regions from calculations
  int i10 = std::fabs(10/lat_res);
  int i171 = std::fabs(171/lat_res);
  double modLevLen = pos_bot-pos_top;
  double invLevLen = 1.0/(2.0*modLevLen);
//Set top/bottom 10 degrees latitude to 0 for PV
    for (int p=0; p<nPlev; p++){
      for (int b=0; b<nLon; b++){
        for (int a=0; a<i10; a++){
          PVMat[p][a][b] = 0.0;
        }
        for (int a=i171; a<nLat; a++){
          PVMat[p][a][b] = 0.0;
        }
      }
    }

  //Calculate integration parts
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        PVtop = PVMat[pos_top][a][b];
        PVbot = PVMat[pos_bot][a][b];
        PVmid = 0.0;
        for (int p=(pos_top+1); p<pos_bot; p++){
          PVmid+=2.0*PVMat[p][a][b];
        }
        IPVMat[a][b] = (PVtop+PVbot+PVmid)*invLevLen;
      }
    }
}


////////////////////////////////////////////////////
//    SECTION: VARIABLE ANOMALY CALCULATIONS      //
////////////////////////////////////////////////////

//Takes vertically averaged variable (PV or Z) and long term 
//daily average and outputs 3 variables: instantaneous 
//anomalies, anomalies with 2-day smoothing, and a normalized
//anomaly (all values below threshold or wrong sign set to 0)
void calcDevs( bool latNorm,
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
              NcVar *lat){

  int nTime,nLat,nLon,avgDay;
  double pi = std::atan(1.)*4.;

  nTime = inIPV->get_dim(0)->size();
  if (is4D=="T"){
    nLat = inIPV->get_dim(2)->size();
    nLon = inIPV->get_dim(3)->size();
  }else{
    nLat = inIPV->get_dim(1)->size();
    nLon = inIPV->get_dim(2)->size();
  }

  DataVector<double> latVec(nLat);
  lat->set_cur((long) 0);
  lat->get(&(latVec[0]),nLat);

  //Vector of instantaneous time axis
  DataVector<double> timeVec(nTime);
  inTime->set_cur((long) 0);
  inTime->get(&(timeVec[0]),nTime);

//Deal with skipped days          
  int d=0;
//  DataVector<double> newTime(nOutTime);

  int leapYear=0;
  int leapMonth=0;
  int leapDay=0;
  int leapHour=0;
  int currAvgIndex=0;
//input instantaneous and average data 
  DataMatrix<double> IPVMat(nLat,nLon);
  DataMatrix<double> avgMat(nLat,nLon);
  DataMatrix<double> devMat(nLat,nLon);
  double num = std::sin(45*pi/180);
  double denom, sineRatio, inputVal;
  for (int t=0; t<nTime; t++){
    if (is4D=="T"){
      inIPV->set_cur(t,pIndex,0,0);
      inIPV->get(&(IPVMat[0][0]),1,1,nLat,nLon);  
    }else{
      inIPV->set_cur(t,0,0);
      inIPV->get(&(IPVMat[0][0]),1,nLat,nLon);
    }
    //check if this time step is a leap day
    ParseTimeDouble(strTimeUnits, strCalendar, timeVec[t], leapYear,\
      leapMonth, leapDay, leapHour);
    if (leapMonth==2 && leapDay == 29){
      std::cout<<"Leap day! Skipping time step."<<std::endl;
    }
    else{
      currAvgIndex = DayInYear(leapMonth,leapDay)-1;
      avgIPV->set_cur(currAvgIndex,0,0);
      avgIPV->get(&(avgMat[0][0]),1,nLat,nLon);
      for (int a=0; a<nLat; a++){
        //Denominator for Z500 anomaly lat normalization
        if (std::fabs(latVec[a])< 5.){
          denom = std::sin(5.*pi/180.);
        }
        else if (std::fabs(latVec[a])>85.){
          denom = std::sin(85. *pi/180.);
        }
        else{
          denom = std::fabs(std::sin(latVec[a]*pi/180.));
        }
        sineRatio = num/denom;
        for (int b=0; b<nLon; b++){
          inputVal = IPVMat[a][b];
          if (ZtoGH=="T"){
            inputVal /= 9.8;
          }
          devMat[a][b] = inputVal-avgMat[a][b];
          if (latNorm){
            devMat[a][b]*=sineRatio;
          }
        }
      }
      outDev->set_cur(d,0,0);
      outDev->put(&(devMat[0][0]),1,nLat,nLon);
      d++;
    }    
  }


}


void calcNormalizedDevs(bool isPV,
                       NcVar * inDev,
                       NcVar * outPosIntDev,
                       NcVar * lat,
                       NcVar * inTime,
                       std::string strTimeUnits,
                       std::string strCalendar,
                       DataMatrix3D<double>threshMat,
                       double minThresh){

  int nLat,nLon,nOutTime;

  nOutTime = inDev->get_dim(0)->size();
  nLat = inDev->get_dim(1)->size();
  nLon = inDev->get_dim(2)->size();
  //lat values
  DataVector<double> latVec(nLat);
  lat->set_cur((long) 0);
  lat->get(&(latVec[0]),nLat);
  //time values
  DataVector<double> timeVals(nOutTime);
  inTime->set_cur((long) 0);
  inTime->get(&(timeVals[0]),nOutTime);

  DataMatrix<double> aDevMat(nLat,nLon);
  DataMatrix<int> posIntDevs(nLat,nLon);
  double invAnom;
  double divDev,pos,neg;
  int threshIndex = 0;
  int dateYear, dateMonth, dateDay, dateHour;
//  int startAvgIndex = 0;
//  int nPastStart = 0;
//  int dPastStart = 0;
  double threshVal;
  if (isPV){
    for (int t=0; t<nOutTime; t++){
      ParseTimeDouble(strTimeUnits,strCalendar,timeVals[t],dateYear,\
         dateMonth,dateDay,dateHour);

      threshIndex = DayInYear(dateMonth,dateDay)-1;
      inDev->set_cur(t,0,0);
      inDev->get(&(aDevMat[0][0]),1,nLat,nLon);
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          if (std::fabs(latVec[a])<25 || std::fabs(latVec[a])>75){
            posIntDevs[a][b] = 0;
          }
          else{
            if (threshMat[threshIndex][a][b] < minThresh){
              threshVal = minThresh;
            }
            else{
              threshVal = threshMat[threshIndex][a][b];
            }
            invAnom = 1./threshVal;
            divDev = aDevMat[a][b]*invAnom;
            //SH: positive anomalies
            if (latVec[a]<0){
              pos = (divDev+std::fabs(divDev))*0.5;
              posIntDevs[a][b] = int(pos);
            }
            //NH: negative anomalies
            else if (latVec[a]>=0){
              neg = (divDev-std::fabs(divDev))*0.5;
              posIntDevs[a][b] = -int(neg);
            }
          }
        }
      }
      outPosIntDev->set_cur(t,0,0);
      outPosIntDev->put(&(posIntDevs[0][0]),1,nLat,nLon);
    }
  }
  else{
    for (int t=0; t<nOutTime; t++){
      ParseTimeDouble(strTimeUnits,strCalendar,timeVals[t],dateYear,\
         dateMonth,dateDay,dateHour);

      threshIndex = DayInYear(dateMonth,dateDay)-1; 
      inDev->set_cur(t,0,0);
      inDev->get(&(aDevMat[0][0]),1,nLat,nLon);
       for (int a=0; a<nLat; a++){
         for (int b=0; b<nLon; b++){
           if (std::fabs(latVec[a])<25. || std::fabs(latVec[a])>75.){
             posIntDevs[a][b] = 0;
           }
           else{
             if (threshMat[threshIndex][a][b] < minThresh){
               threshVal = minThresh;
             }
             else{
               threshVal = threshMat[threshIndex][a][b];
             }
             invAnom = 1./threshVal;
             pos = aDevMat[a][b]*invAnom;
             posIntDevs[a][b] = (int)((pos + std::fabs(pos))*0.5);
          }
        }
      }
      outPosIntDev->set_cur(t,0,0);
      outPosIntDev->put(&(posIntDevs[0][0]),1,nLat,nLon);
    }
  }

  std::cout<<"Wrote integer values to file."<<std::endl;
}

/*
void stdDev(DataMatrix3D<double>inDevs,
              int nTime,
              int nLat,
              int nLon,
              DataMatrix<double> & outStdDev){
  //average the deviations along the time axis
  double sigSum;
  double dev;
  double variance;
  
  for (int a=0; a<nLat; a++){
    for (int b=0; b<nLon; b++){
      sigSum = 0;
      for (int t=0; t<nTime; t++){
        dev = inDevs[t][a][b];
        sigSum += (dev*dev);
      }
      variance = sigSum/nTime;
      outStdDev[a][b] = std::sqrt(variance);
    //  if (a==32 && b== 50){
    //    std::cout<<"Standard deviation value: "<<std::sqrt(variance)<<std::endl;
   //   }
    }
  }
}
*/
/*
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
              DataMatrix3D threshmat){

  int nTime,nLat,nLon,nSteps,avgDay,nOutTime;
  double tRes;
  double pi = std::atan(1)*4.;

  nTime = inGH->get_dim(0)->size();
  nLat = inGH->get_dim(1)->size();
  nLon = inGH->get_dim(2)->size();

//input GH
  DataMatrix3D<double> GHMat(nTime,nLat,nLon);
  inGH->set_cur(0,0,0);
  inGH->get(&(GHMat[0][0][0]),nTime,nLat,nLon);

  DataVector<double> timeVec(nTime);
  inTime->set_cur((long) 0);
  inTime->get(&(timeVec[0]),nTime);

  std::string strTimeUnits = inTime->get_att("units")->as_string(0);
  std::string strCalendar = inTime->get_att("calendar")->as_string(0);

  if ((strTimeUnits.length() >= 11) && \
    (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)){
    tRes = timeVec[1]-timeVec[0];
  }
  else{
    tRes = (timeVec[1]-timeVec[0])/24.0;
  }
  nSteps = 1/tRes;

//avg GH
  avgDay = avgGH->get_dim(0)->size();
  DataMatrix3D<double> avgMat(avgDay,nLat,nLon);
  avgGH->set_cur(0,0,0);
  avgGH->get(&(avgMat[0][0][0]),avgDay,nLat,nLon);

  DataVector<int> avgDayVec(avgDay);
  avgTime->set_cur((long) 0);
  avgTime->get(&(avgDayVec[0]),avgDay);

//Latitude values 

  DataVector<double> latVec(nLat);
  lat->set_cur((long) 0);
  lat->get(&(latVec[0]),nLat);
//Matrix for output data
//Eliminate one day if contains Feb 29
  if (leap){
    nOutTime = nTime-nSteps;
  }
  else{
    nOutTime = nTime;
  }

  DataMatrix3D<double> devMat(nOutTime,nLat,nLon);
  DataMatrix3D<double> aDevMat(nOutTime,nLat,nLon);

//Number of days in IPV
  int nDays = nTime*tRes;
  //std::cout<<"There are "<<nDays<<" days in file."<<std::endl;


//Deal with skipped days          
  int d=0;
  DataVector<double> newTime(nOutTime);

  int leapYear=0;
  int leapMonth=0;
  int leapDay=0;
  int leapHour=0;

  //std::cout<<"Starting avg index is "<<startAvgIndex<<std::endl;

  for (int t=0; t<nTime; t++){
    if (leap){
      ParseTimeDouble(strTimeUnits, strCalendar, timeVec[t], leapYear,\
        leapMonth, leapDay, leapHour);
      if (leapMonth==2 && leapDay == 29){
        //std::cout<<"Leap day! Skipping day."<<std::endl;
        t+=nSteps;
      }
    }
    if (t>=nTime){
      break;
    }
    int nDayIncrease = d/nSteps;
   // std::cout<<"Number of days increased since start is "<<nDayIncrease<<std::endl;
    int currAvgIndex = startAvgIndex + nDayIncrease;
    //std::cout<<"Avg index:"<<currAvgIndex<<std::endl;
    if (currAvgIndex>364){
      currAvgIndex-=365;
    //  std::cout<<"Going back to beginning of average index."<<std::endl;
    }
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        devMat[t][a][b] = GHMat[t][a][b]-avgMat[currAvgIndex][a][b];
      }
    }
    newTime[d] = timeVec[t];
    //std::cout<<"d,t:"<<d<<","<<t<<std::endl;
    d++;
  }
  outTime->set_cur((long) 0);
  outTime->put(&(newTime[0]),nOutTime);

  std::cout<<"About to implement smoothing."<<std::endl;
  double div = (double) 2*nSteps;
  double invDiv = 1.0/div;

  //implement 2-day smoothing
  for (int t=0; t<nOutTime; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        if (t<2*nSteps){
          aDevMat[t][a][b] = devMat[t][a][b];
        }
        else{
          for (int n=0; n<2*nSteps; n++){
            aDevMat[t][a][b]+=devMat[t-n][a][b];
          }
          aDevMat[t][a][b] = aDevMat[t][a][b]*invDiv;
        }
      }
    }
  }

  std::cout<<"Finished smoothing."<<std::endl;

  double num = std::sin(45*pi/180);
  double denom;
  double sineRatio;
  //Multiply values by sine factor
  for (int t=0; t<nOutTime; t++){
    for (int a=0; a<nLat; a++){
      if ((latVec[a] < 30. && latVec[a] > -30.)\
        || latVec[a] > 75. || latVec[a] < -75.){
        sineRatio = 0.;
      }
      else{
        denom = std::fabs(std::sin(latVec[a]*pi/180.));
        sineRatio = num/denom;
      }
      for (int b=0; b<nLon; b++){
        devMat[t][a][b]*=sineRatio;
        aDevMat[t][a][b]*=sineRatio;
      }
    }
  }

  outDev->set_cur(0,0,0);
  outDev->put(&(devMat[0][0][0]),nOutTime,nLat,nLon);
  std::cout<<"Wrote devs to file."<<std::endl;

  outADev->set_cur(0,0,0);
  outADev->put(&(aDevMat[0][0][0]),nOutTime,nLat,nLon);
  std::cout<<"Wrote smoothed devs to file."<<std::endl;
//Divide matrix by GH anomaly value 

  DataMatrix3D<int> posIntDevs(nOutTime,nLat,nLon);
  DataMatrix<double> stdDevs(nLat,nLon);
  stdDev(aDevMat,nOutTime,nLat,nLon,stdDevs);

  double invAnom;
  double pos;
 
  for (int a=0; a<nLat; a++){
    for (int b=0; b<nLon; b++){
      if (std::fabs(latVec[a])<30. || std::fabs(latVec[a])>75.){
        invAnom = 0.;
      }
      else{
       invAnom = 1./GHAnom;;
      }
      for (int t=0; t<nOutTime; t++){
        pos = aDevMat[t][a][b]*invAnom;
        posIntDevs[t][a][b] =(int)((pos + std::fabs(pos))*0.5);
        }
      }
    }
 
  stdDevVar->set_cur(0,0);
  stdDevVar->put(&(stdDevs[0][0]),nLat,nLon);
 
  outIntDev->set_cur(0,0,0);
  outIntDev->put(&(posIntDevs[0][0][0]),nOutTime,nLat,nLon);
  std::cout<<"Wrote integer values to file."<<std::endl;
}



*/



void PV_calc2(
        int nPlev,
        int nLat,
        int nLon,
	DataMatrix3D<double>UMat,
	DataMatrix3D<double>VMat,
	DataMatrix3D<double> PTMat,
	DataMatrix3D<double> RVMat,
        DataVector<double>pVec,	
	DataVector<double> coriolis,
        DataVector<double>cosphi,
	double dphi,
	double dlambda,
        double lat_res,
        double lon_res,
	DataMatrix3D<double> &PVMat,
  DataMatrix3D<double> &dpt_dp,
  DataMatrix3D<double> &du_dp,
  DataMatrix3D<double> &dv_dp,
  DataMatrix3D<double> &dpt_dphi,
  DataMatrix3D<double> &dpt_dl 
){
  double invdp,invdp1,invdp2;
  double invdphi= 1.0/(2.0*dphi);
  double invdlambda = 1.0/(2.0*dlambda);
  double radius = 6371000.0;
  double coef1,coef2,corvar;
  double U1,U2,U3,V1,V2,V3;


  //Matrices for the partials
  //PT, U, V WRT P

  invdp1 = 1.0/(2.0*(pVec[1]-pVec[0]));
  invdp2 = 1.0/(2.0*(pVec[nPlev-1]-pVec[nPlev-2]));

    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        U1=UMat[2][a][b];
        if (std::fabs(U1)>1000.0){
          U1 = replaceMissingFloat(a,b,2,1000.0,UMat,nLat,nLon);
        }
        U2=UMat[1][a][b];
        if (std::fabs(U2)>1000.0){
          U2 = replaceMissingFloat(a,b,1,1000.0,UMat,nLat,nLon);
        }
        U3=UMat[0][a][b];
        if (std::fabs(U3)>1000.0){
          U3 = replaceMissingFloat(a,b,0,1000.0,UMat,nLat,nLon);
        }
       
        V1=VMat[2][a][b];
        if (std::fabs(V1)>1000.0){
          V1 = replaceMissingFloat(a,b,2,1000.0,VMat,nLat,nLon);
        }
        V2=VMat[1][a][b];
        if (std::fabs(V2)>1000.0){
          V2 = replaceMissingFloat(a,b,1,1000.0,VMat,nLat,nLon);
        }
        V3=VMat[0][a][b];
        if (std::fabs(V3)>1000.0){
          V3 = replaceMissingFloat(a,b,0,1000.0,VMat,nLat,nLon);
        }
      //0 case
        dpt_dp[0][a][b] = (-PTMat[2][a][b]+4.0*PTMat[1][a][b]-3.0*PTMat[0][a][b])*invdp1;
        du_dp[0][a][b] = (-U1+4.0*U2-3.0*U3)*invdp1;
        dv_dp[0][a][b] = (-V1+4.0*V2-3.0*V3)*invdp1;          
      //end case
        U1=UMat[nPlev-1][a][b];
        if (std::fabs(U1)>1000.0){
          U1 = replaceMissingFloat(a,b,(nPlev-1),1000.0,UMat,nLat,nLon);
        }
        U2=UMat[nPlev-2][a][b];
        if (std::fabs(U2)>1000.0){
          U2 = replaceMissingFloat(a,b,(nPlev-2),1000.0,UMat,nLat,nLon);
        }
        U3=UMat[nPlev-3][a][b];
        if (std::fabs(U3)>1000.0){
          U3 = replaceMissingFloat(a,b,(nPlev-3),1000.0,UMat,nLat,nLon);
        }
       
        V1=VMat[nPlev-1][a][b];
        if (std::fabs(V1)>1000.0){
          V1 = replaceMissingFloat(a,b,(nPlev-1),1000.0,VMat,nLat,nLon);
        }
        V2=VMat[nPlev-2][a][b];
        if (std::fabs(V2)>1000.0){
          V2 = replaceMissingFloat(a,b,nPlev-2,1000.0,VMat,nLat,nLon);
        }
        V3=VMat[nPlev-3][a][b];
        if (std::fabs(V3)>1000.0){
          V3 = replaceMissingFloat(a,b,nPlev-3,1000.0,VMat,nLat,nLon);
        }
        dpt_dp[nPlev-1][a][b] = (3.0*PTMat[nPlev-1][a][b]-4.0*PTMat[nPlev-2][a][b]\
          +PTMat[nPlev-3][a][b])*invdp2;
        du_dp[nPlev-1][a][b] = (3.0*U1-4.0*U2+U3)*invdp2;
        dv_dp[nPlev-1][a][b] = (3.0*V2-4.0*V2+V3)*invdp2;
        for (int p=1; p<(nPlev-1); p++){
          U1=UMat[p+1][a][b];
          if (std::fabs(U1)>1000.0){
            U1 = replaceMissingFloat(a,b,(p+1),1000.0,UMat,nLat,nLon);
          }
          U2=UMat[p-1][a][b];
          if (std::fabs(U2)>1000.0){
            U2 = replaceMissingFloat(a,b,(p-1),1000.0,UMat,nLat,nLon);
          }
       
          V1=VMat[p+1][a][b];
          if (std::fabs(V1)>1000.0){
            V1 = replaceMissingFloat(a,b,(p-1),1000.0,VMat,nLat,nLon);
          }
          V2=VMat[p-1][a][b];
          if (std::fabs(V2)>1000.0){
            V2 = replaceMissingFloat(a,b,(p-1),1000.0,VMat,nLat,nLon);
          }
          invdp = 1.0/(2.0*(pVec[p+1]-pVec[p]));
          dpt_dp[p][a][b] = (PTMat[p+1][a][b]-PTMat[p-1][a][b])*invdp;
          du_dp[p][a][b] = (U1-U2)*invdp;
          dv_dp[p][a][b] = (V1-V2)*invdp;
        }
      }
    }
  
  //PT WRT PHI
  //end cases
    for (int p=0; p<nPlev; p++){
      for (int b=0; b<nLon; b++){
        dpt_dphi[p][0][b]=(-PTMat[p][2][b]+4.0*PTMat[p][1][b]\
          -3.0*PTMat[p][0][b])*invdphi;
        dpt_dphi[p][nLat-1][b]=(3.0*PTMat[p][nLat-1][b]-4.0*PTMat[p][nLat-2][b]\
          +PTMat[p][nLat-2][b])*invdphi;
        for (int a=1; a<(nLat-1); a++){
          dpt_dphi[p][a][b]=(PTMat[p][a+1][b]-PTMat[p][a-1][b])*invdphi;
        }
      }
    }

  //PT WRT LAMBDA
  //end cases
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        dpt_dl[p][a][0]=(PTMat[p][a][1]-PTMat[p][a][nLon-1])*invdlambda;
        dpt_dl[p][a][nLon-1]=(PTMat[p][a][nLon-2]-PTMat[p][a][0])*invdlambda;
        for (int b=1; b<(nLon-1); b++){
          dpt_dl[p][a][b]=(PTMat[p][a][b+1]-PTMat[p][a][b-1])*invdlambda;
        }
      }
    }
  coef2 = 1.0/radius;
  //PV Calculation!
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        coef1 = 1.0/(radius*cosphi[a]);
        corvar=coriolis[a];
        for (int b=0; b<nLon; b++){
          PVMat[p][a][b] = 9.80616*(coef1*dv_dp[p][a][b]*dpt_dl[p][a][b]\
           -coef2*du_dp[p][a][b]*dpt_dphi[p][a][b]\
           -(corvar+RVMat[p][a][b])*dpt_dp[p][a][b]);
           if (std::fabs(PVMat[p][a][b])>10){
             std::cout<<"Questionable value at p="<<p<<", a="<<a<<", b="<<b<<std::endl;
           }
        }
      }
    }
} 

