/////////////////////////////////////////////
///
///          \file CLIVAR_block_utilities.cpp
///
///          \author Marielle Pinheiro
///          
///          \version March 1, 2015
///

#include "CLIVAR_block_utilities.h"
#include "NetCDFUtilities.h"
#include "netcdfcpp.h"
#include "DataVector.h"
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


//Function to get number of day in year

int DayInYear(int nMonth, int nDay){
  int day=0;
  std::cout<<"nMonth is "<<nMonth<<std::endl;
  if (nMonth>1){
    for (int x=1; x<nMonth; x++){
      std::cout<<"Month is currently "<<x;
      if (x==2){
        day += 28;
        std::cout<<". Adding 28 days.";
      }
      else if (x==4 || x==6 || x==9 || x==11){
        day +=30;
        std::cout<<". Adding 30 days.";
      }
      else{
        day +=31;
        std::cout<<". Adding 31 days.";
      }
      std::cout<<" Total number of days currently "<<day<<std::endl;
    }
  }
  day+=nDay;
  std::cout<<"Final day value is "<<day<<std::endl;
  return(day);
}



//Copied from DetectCyclones
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

	} else {
		_EXCEPTION1("Unknown calendar type \"%s\"", strTimeCalendar.c_str());
	}
/*
	Time time(Time::CalendarStandard);
	time.FromFormattedString("1800-01-01 00:00:00");
	printf("%1.15e %i\n", 3600.0 * 1577832.0, (int)(3600.0 * 1577832.0));
	time.AddHours(1577832);

	Announce("Time (YMDS): %i %i %i %i",
			time.GetYear(),
			time.GetMonth(),
			time.GetDay(),
			time.GetSecond());

	_EXCEPTION();
*/
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

		Announce("Time (YMDS): %i %i %i %i",
				time.GetYear(),
				time.GetMonth(),
				time.GetDay(),
				time.GetSecond());


		nDateYear = time.GetYear();
		nDateMonth = time.GetMonth();
		nDateDay = time.GetDay();
		nDateHour = time.GetSecond() / 3600;

		//printf("%s\n", strSubStr.c_str());

	// Time format is "hours since ..."
	} else if (
	    (strTimeUnits.length() >= 12) &&
	    (strncmp(strTimeUnits.c_str(), "hours since ", 12) == 0)
	) {
		std::string strSubStr = strTimeUnits.substr(12);
		Time time(cal);
		time.FromFormattedString(strSubStr);

		time.AddHours(static_cast<int>(dTime));

		Announce("Time (YMDS): %i %i %i %i",
				time.GetYear(),
				time.GetMonth(),
				time.GetDay(),
				time.GetSecond());

		nDateYear = time.GetYear();
		nDateMonth = time.GetMonth();
		nDateDay = time.GetDay();
		nDateHour = time.GetSecond() / 3600;

		//printf("%s\n", strSubStr.c_str());

	} else {
		_EXCEPTIONT("Unknown \"time::units\" format");
	}
	//_EXCEPTION();
}






///////////////////////////////////////////////////////////////////////////////////////
//Function to interpolate variables from hybrid levels to pressure levels
void interpolate_lev(NcVar *var, 
                     NcVar *hyam, 
                     NcVar *hybm, 
                     NcVar *ps, 
                     NcVar *pLev,
                     NcVar *NewVar 
){
  //Array dimensions
  int nTime = var->get_dim(0)->size();
  int nLev = var->get_dim(1)->size();
  int nLat = var->get_dim(2)->size();
  int nLon = var->get_dim(3)->size();
  int npLev = pLev->get_dim(0)->size();

 // std::cout<<"Dimensions: Time: "<<nTime<<", hybrid level: "<<nLev\
    << " lat: " << nLat << " lon: " << nLon << " plev: " << npLev<<std::endl;

  //Matrix to store PS
  DataMatrix3D <double> matPS(nTime, nLat, nLon);
  ps->set_cur(0, 0, 0);
  ps->get(&(matPS[0][0][0]), nTime, nLat, nLon);

  //Matrix to store input variable data
  DataMatrix4D<double> matVar(nTime, nLev, nLat, nLon);
  var->set_cur(0, 0, 0, 0);
  var->get(&(matVar[0][0][0][0]), nTime, nLev, nLat, nLon);
  
  //Matrix to store output variable data
  DataMatrix4D<double> matVarOut(nTime, npLev, nLat, nLon);

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
  
  //Loop over input data and interpolate to output var
  //Need to deal with areas of 0!
  for (int t=0; t<nTime; t++){
    for (int l=0; l<(nLev-1); l++){
      double A1 = vecHyam[l];
      double B1 = vecHybm[l];
      double A2 = vecHyam[l+1];
      double B2 = vecHybm[l+1];
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
	  int interp_check = 0;
          double p1 = 100000.0 * A1 + matPS[t][a][b] * B1;
          double p2 = 100000.0 * A2 + matPS[t][a][b] * B2;
          for (int p=0; p<npLev; p++){
            if (p1<vecpLev[p] && p2>=vecpLev[p]){
	      interp_check = 1;
              double weight = ((vecpLev[p]-p1)/(p2-p1));
              matVarOut[t][p][a][b] = weight*matVar[t][l+1][a][b]
                               + (1.0-weight)*matVar[t][l][a][b];
            }
          }
        }
      }
    }
  }
  std::cout<<"Finished interpolating variable.\n";
  NewVar->set_cur(0, 0, 0, 0);
  NewVar->put(&(matVarOut[0][0][0][0]), nTime, npLev, nLat, nLon);
  CopyNcVarAttributes(var, NewVar);
}

////////////////////////////////////////////////////////////////////////

//Function to copy dimension variables to outfile
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

///////////////////////////////////////////////////////////////////////

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
){ 
  double radius = 6371000.0;
  double pi = 4.0*std::atan(1.0);
  double radian = 180.0/pi;
  double sigma = 7.2921*std::pow(10.0, -5.0);
//  std::cout<<"Check: radius is "<<radius<<", pi is "<<pi<<", radian is "\
    <<radian<<", and sigma is "<<sigma<<std::endl;
  //std::cout<<"Initiating var calcs."<<std::endl;
  int nLat = lat->get_dim(0)->size();
  DataVector<double> latVec(nLat);
  lat->set_cur((long) 0);
  lat->get(&(latVec[0]), nLat);

  int nLon = lon->get_dim(0)->size();
  DataVector<double> lonVec(nLon);
  lon->set_cur((long) 0);
  lon->get(&(lonVec[0]), nLon);

  int np = plev->get_dim(0)->size();
  DataVector<double> pVec(np);
  plev->set_cur((long) 0);
  plev->get(&(pVec[0]), np);
  
 
  lat_res = latVec[2]-latVec[1];
  lon_res = lonVec[2]-lonVec[1];
  dphi = lat_res/radian;
  dlambda = lon_res/radian;
  p_res = pVec[2]-pVec[1];
  for (int i=0; i<nLat; i++){
    double sinphi = std::sin(latVec[i]/radian);
    coriolis[i] = 2.0 * sigma * sinphi;
    cosphi[i] = std::cos(latVec[i]/radian);
  }
}
/////////////////////////////////////////////////////////////////

//Function that calculates PT
void PT_calc(
	NcVar *T, 
	NcVar *pLev, 
	DataMatrix4D<double> &PTMat
){
  int nTime = T->get_dim(0)->size();
  int nPlev = T->get_dim(1)->size();
  int nLat = T->get_dim(2)->size();
  int nLon = T->get_dim(3)->size();

  //INPUTS
  DataMatrix4D<double> TMat(nTime, nPlev, nLat, nLon);
  DataVector<double> pVec(nPlev);

  T->set_cur(0,0,0,0);
  T->get(&(TMat[0][0][0][0]), nTime, nPlev, nLat, nLon);

  pLev->set_cur((long) 0);
  pLev->get(&(pVec[0]), nPlev);

  //OUTPUT: PT
//  DataMatrix4D<double> PTMat(nTime, nPlev, nLat, nLon);
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          double p_frac = 100000.0/pVec[p];
          double exp = 287.0/1004.5;
          PTMat[t][p][a][b] = TMat[t][p][a][b]*std::pow(p_frac, exp);
        }
      }
    }
  }
 // PT->set_cur(0,0,0,0);
 // PT->put(&(PTMat[0][0][0][0]), nTime, nPlev, nLat, nLon);
  std::cout<<"Finished calculating PT."<<std::endl;
}
///////////////////////////////////////////////////////////////////////

//Function that calculates relative vorticity
void rVort_calc(
	NcVar *U,
	NcVar *V,
	double dphi,
	double dlambda,
	DataVector<double> cosphi,
	DataMatrix4D<double> & RVMat
){
  int nTime = U->get_dim(0)->size();
  int nPlev = U->get_dim(1)->size();
  int nLat = U->get_dim(2)->size();
  int nLon = U->get_dim(3)->size();

 double radius = 6371000.0;
//INPUTS
  DataMatrix4D<double> VMat(nTime, nPlev, nLat, nLon);
  V->set_cur(0,0,0,0);
  V->get(&(VMat[0][0][0][0]), nTime, nPlev, nLat, nLon);

  DataMatrix4D<double> UMat(nTime, nPlev, nLat, nLon);
  U->set_cur(0,0,0,0);
  U->get(&(UMat[0][0][0][0]),nTime, nPlev, nLat, nLon);
//OUTPUT: Partial derivatives
//U WRT PHI
  DataMatrix4D<double> dUdphi(nTime, nPlev, nLat, nLon);
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int b=0; b<nLon; b++){
        dUdphi[t][p][0][b]=(-UMat[t][p][2][b]*cosphi[2]\
          +4.0*UMat[t][p][1][b]*cosphi[1]\
          -3.0*UMat[t][p][0][b]*cosphi[0])/(2.0*dphi);
        dUdphi[t][p][nLat-1][b]=(3.0*UMat[t][p][nLat-1][b]*cosphi[nLat-1]\
          -4.0*UMat[t][p][nLat-2][b]*cosphi[nLat-2]\
          +UMat[t][p][nLat-3][b]*cosphi[nLat-3])/(2.0*dphi);
        for (int a=1; a<(nLat-1); a++){
          dUdphi[t][p][a][b]=(UMat[t][p][a+1][b]*cosphi[a+1]\
           -UMat[t][p][a-1][b]*cosphi[a-1])/(2.0*dphi);
        }
      }
    }
  }
  /*for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int a=1; a<(nLat-1); a++){
        for (int b=0; b<nLon; b++){
          dUdphi[t][p][a][b]=(UMat[t][p][a+1][b]*cosphi[a+1]\
           -UMat[t][p][a-1][b]*cosphi[a-1])/(2.0*dphi);
        }
      }
    }
  }*/


 /* //Lat end cases: set to 0 because of pole singularities causing error
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int b=0; b<nLon; b++){
        dUdphi[t][p][0][b] = 0.0;
        dUdphi[t][p][nLat-1][b] = 0.0;
      }
    }
  }*/

  //V WRT LAMBDA
  DataMatrix4D<double> dVdl(nTime, nPlev, nLat, nLon);
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        dVdl[t][p][a][0]=(VMat[t][p][a][1]-VMat[t][p][a][nLon-1])/(2.0*dlambda);
        dVdl[t][p][a][nLon-1]=(VMat[t][p][a][0]-VMat[t][p][a][nLon-2])/(2.0*dlambda);
        for (int b=1; b<(nLon-1); b++){
          dVdl[t][p][a][b]=(VMat[t][p][a][b+1]-VMat[t][p][a][b-1])/(2.0*dlambda);
        }
      }
    }
  }
  //DataMatrix4D<double>RVMat(nTime, nPlev, nLat, nLon);
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          RVMat[t][p][a][b] = (1.0/(radius*cosphi[a]))*(dVdl[t][p][a][b]-dUdphi[t][p][a][b]);
        }
      }
    }
  }

  //Lat end cases: set to 0 because of pole singularities causing error
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int b=0; b<nLon; b++){
        RVMat[t][p][0][b] = 0.0;
        RVMat[t][p][nLat-1][b] = 0.0;
      }
    }
  }

  std::cout<<"Finished calculating relative vorticity."<<std::endl;
 // vorticity->set_cur(0,0,0,0);
 // vorticity->put(&(RVMat[0][0][0][0]),nTime, nPlev, nLat, nLon);
}

///////////////////////////////////////////////////////////////////////////////////////////////

void PV_calc(
	NcVar *U, 
	NcVar *V, 
	DataMatrix4D<double> PTMat,
	DataMatrix4D<double> RVMat,
        NcVar *pVals, 	
	DataVector<double> coriolis,
        DataVector<double>cosphi,
	double dphi,
	double dlambda,
        double lat_res,
        double lon_res,
	NcVar *PV,
        NcVar *intPV){

  int nTime = U->get_dim(0)->size();
  int nPlev = U->get_dim(1)->size();
  int nLat = U->get_dim(2)->size();
  int nLon = U->get_dim(3)->size();

  //Input matrices
  DataMatrix4D<double> UMat(nTime, nPlev, nLat, nLon);
  DataMatrix4D<double> VMat(nTime, nPlev, nLat, nLon);
//  DataMatrix4D<double> PTMat(nTime, nPlev, nLat, nLon);
//  DataMatrix4D<double> RVMat(nTime, nPlev, nLat, nLon);

  //Load data
  U->set_cur(0,0,0,0); 
  U->get(&(UMat[0][0][0][0]),nTime, nPlev, nLat, nLon);

  V->set_cur(0,0,0,0);
  V->get(&(VMat[0][0][0][0]), nTime, nPlev, nLat, nLon);

//  PT->set_cur(0,0,0,0);
//  PT->get(&(PTMat[0][0][0][0]), nTime, nPlev, nLat, nLon);

//  rVort->set_cur(0,0,0,0);
//  rVort->get(&(RVMat[0][0][0][0]), nTime, nPlev, nLat, nLon);

  //Pressure axis values
  DataVector<double> pVec(nPlev);
  pVals->set_cur((long) 0);
  pVals->get(&(pVec[0]), nPlev);

  //Matrices for the partials
  //PT, U, V WRT P
  DataMatrix4D<double> dpt_dp(nTime, nPlev, nLat, nLon);
  DataMatrix4D<double> du_dp(nTime, nPlev, nLat, nLon);
  DataMatrix4D<double> dv_dp(nTime, nPlev, nLat, nLon);

  double dp1 = std::fabs(pVec[1]-pVec[0]);
  double dp2 = std::fabs(pVec[nPlev-1]-pVec[nPlev-2]);
  for (int t=0; t<nTime; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
      //0 case
        dpt_dp[t][0][a][b] = (-PTMat[t][2][a][b]+4.0*PTMat[t][1][a][b]-3.0*PTMat[t][0][a][b])/(2.0*dp1);
        du_dp[t][0][a][b] = (-UMat[t][2][a][b]+4.0*UMat[t][1][a][b]-3.0*UMat[t][0][a][b])/(2.0*dp1);
        dv_dp[t][0][a][b] = (-VMat[t][2][a][b]+4.0*VMat[t][1][a][b]-3.0*VMat[t][0][a][b])/(2.0*dp1);           //end case
        dpt_dp[t][nPlev-1][a][b] = (3.0*PTMat[t][nPlev-1][a][b]-4.0*PTMat[t][nPlev-2][a][b]\
          +PTMat[t][nPlev-3][a][b])/(2.0*dp2);
        du_dp[t][nPlev-1][a][b] = (3.0*UMat[t][nPlev-1][a][b]-4.0*UMat[t][nPlev-2][a][b]\
          +UMat[t][nPlev-3][a][b])/(2.0*dp2);
        dv_dp[t][nPlev-1][a][b] = (3.0*VMat[t][nPlev-1][a][b]-4.0*VMat[t][nPlev-2][a][b]\
          +VMat[t][nPlev-3][a][b])/(2.0*dp2);
        for (int p=1; p<(nPlev-1); p++){
          double dp = std::fabs(pVec[p+1]-pVec[p]);
          dpt_dp[t][p][a][b] = (PTMat[t][p+1][a][b]-PTMat[t][p-1][a][b])/(2.0*dp);
          du_dp[t][p][a][b] = (UMat[t][p+1][a][b]-UMat[t][p-1][a][b])/(2.0*dp);
          dv_dp[t][p][a][b] = (VMat[t][p+1][a][b]-VMat[t][p-1][a][b])/(2.0*dp);
        }
      }
    }
  }
  
  //PT WRT PHI
  DataMatrix4D<double> dpt_dphi(nTime, nPlev, nLat, nLon);
  //end cases
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int b=0; b<nLon; b++){
        dpt_dphi[t][p][0][b]=(-PTMat[t][p][2][b]+4.0*PTMat[t][p][1][b]\
          -3.0*PTMat[t][p][0][b])/(2.0*dphi);
        dpt_dphi[t][p][nLat-1][b]=(3.0*PTMat[t][p][nLat-1][b]-4.0*PTMat[t][p][nLat-2][b]\
          +PTMat[t][p][nLat-2][b])/(2.0*dphi);
        for (int a=1; a<(nLat-1); a++){
          dpt_dphi[t][p][a][b]=(PTMat[t][p][a+1][b]-PTMat[t][p][a-1][b])/(2.0*dphi);
        }
      }
    }
  }

  //PT WRT LAMBDA
  DataMatrix4D<double> dpt_dl(nTime, nPlev, nLat, nLon);
  //end cases
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        dpt_dl[t][p][a][0]=(PTMat[t][p][a][1]-PTMat[t][p][a][nLon-1])/(2.0*dlambda);
        dpt_dl[t][p][a][nLon-1]=(PTMat[t][p][a][nLon-2]-PTMat[t][p][a][0])/(2.0*dlambda);
        for (int b=1; b<(nLon-1); b++){
          dpt_dl[t][p][a][b]=(PTMat[t][p][a][b+1]-PTMat[t][p][a][b-1])/(2.0*dlambda);
        }
      }
    }
  }
 
  DataMatrix4D<double> PVMat(nTime, nPlev, nLat, nLon);
  //PV Calculation!
  double radius = 6371000.0;
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        double cosvar=cosphi[a];
        double corvar=coriolis[a];
        for (int b=0; b<nLon; b++){
          PVMat[t][p][a][b] = 9.80616*((1.0/(radius*cosvar))*dv_dp[t][p][a][b]*dpt_dl[t][p][a][b]\
           -(1.0/radius)*du_dp[t][p][a][b]*dpt_dphi[t][p][a][b]\
           -(corvar+RVMat[t][p][a][b])*dpt_dp[t][p][a][b]);
        }
      }
    }
  }
  std::cout<<"Calculated potential vorticity."<<std::endl;           
  PV->set_cur(0,0,0,0);
  PV->put(&(PVMat[0][0][0][0]),nTime,nPlev,nLat,nLon);

  //Integrate PV over upper troposphere
  int pos_top;
  int pos_bot;

  for (int x=0; x<nPlev; x++){
    std::cout<<"pvec value is currently "<<pVec[x]<<std::endl;
    if (std::fabs(pVec[x]-15000.0)<0.0001){
      pos_top = x;
      std::cout<<"diff is "<<std::fabs(pVec[x]-15000.0)<<", 150 mb position is at x="<<x<<std::endl;
    }
    if (abs(pVec[x]-50000.0)<0.0001){
      pos_bot = x;
      std::cout<<"500 mb position is at x="<<x<<std::endl;
    }
  }
 //if level axis starts at ground, reverse positions
  if (pos_top>pos_bot){
    std::cout<<"Reversing axis positions."<<std::endl;
    int temp = pos_bot;
    pos_bot = pos_top;
    pos_top = temp;
  }


  DataMatrix3D<double> iPVMat(nTime, nLat, nLon);

  double PVbot;
  double PVmid;
  double PVtop;
  
//Eliminate polar regions from calculations
  int i10 = std::fabs(10/lat_res);
  int i171 = std::fabs(171/lat_res);
  double modLevLen = pos_bot-pos_top;
//  std::cout<<"Length of modified pressure axis is "<<int_length<<std::endl;

//Set top/bottom 10 degrees latitude to 0 for PV
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int b=0; b<nLon; b++){
        for (int a=0; a<i10; a++){
          PVMat[t][p][a][b] = 0.0;
        }
        for (int a=i171; a<nLat; a++){
          PVMat[t][p][a][b] = 0.0;
        }
      }
    }
  }

  //Calculate integration parts
  for (int t=0; t<nTime; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        PVtop = PVMat[t][pos_top][a][b];
        PVbot = PVMat[t][pos_bot][a][b];
        PVmid = 0.0;
        for (int p=(pos_top+1); p<pos_bot; p++){
          PVmid+=2.0*PVMat[t][p][a][b];
        }
        iPVMat[t][a][b] = (PVtop+PVbot+PVmid)/(2.0*modLevLen);
      }
    }
  }
 
  intPV->set_cur(0,0,0,0);
  intPV->put(&(iPVMat[0][0][0]),nTime,nLat,nLon);
  std::cout<<"Finished integrating PV."<<std::endl;
} 
//////////////////////////////////////////////////////////////////////////////////////////////////


