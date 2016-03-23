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

	} else if (
                (strTimeCalendar.length() >=9) &&
                (strncmp(strTimeCalendar.c_str(), "gregorian",9)==0)
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

	/*	Announce("Time (YMDS): %i %i %i %i",
				time.GetYear(),
				time.GetMonth(),
				time.GetDay(),
				time.GetSecond());*/


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
              //  printf("Debug: dTime is %10f \n",dTime);
		time.AddHours(static_cast<int>(dTime));

	/*	Announce("Time (YMDS): %i %i %i %i",
				time.GetYear(),
				time.GetMonth(),
				time.GetDay(),
				time.GetSecond());*/

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

  int nTime,nLev,nLat,nLon,npLev;
  double A1,A2,B1,B2,p1,p2,weight;

  //Array dimensions
  nTime = var->get_dim(0)->size();
  nLev = var->get_dim(1)->size();
  nLat = var->get_dim(2)->size();
  nLon = var->get_dim(3)->size();
  npLev = pLev->get_dim(0)->size();

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
  for (int t=0; t<nTime; t++){
    for (int l=0; l<(nLev-1); l++){
      A1 = vecHyam[l];
      B1 = vecHybm[l];
      A2 = vecHyam[l+1];
      B2 = vecHybm[l+1];
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
	  int interp_check = 0;
          p1 = 100000.0 * A1 + matPS[t][a][b] * B1;
          p2 = 100000.0 * A2 + matPS[t][a][b] * B2;
          for (int p=0; p<npLev; p++){
            if (p1<vecpLev[p] && p2>=vecpLev[p]){
	      interp_check = 1;
              weight = ((vecpLev[p]-p1)/(p2-p1));
              matVarOut[t][p][a][b] = weight*matVar[t][l+1][a][b]
                               + (1.0-weight)*matVar[t][l][a][b];
            }
          }
        }
      }
    }
  }
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
/////////////////////////////////////////////////////////////////

//Function that calculates PT
void PT_calc(
	NcVar *T, 
	NcVar *pLev, 
	DataMatrix4D<double> &PTMat
){

  int nTime,nPlev,nLat,nLon;
  double pFrac;
  double exp = 287.0/1004.5;

  nTime = T->get_dim(0)->size();
  nPlev = T->get_dim(1)->size();
  nLat = T->get_dim(2)->size();
  nLon = T->get_dim(3)->size();

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
      pFrac = 100000.0/pVec[p];
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          PTMat[t][p][a][b] = TMat[t][p][a][b]*std::pow(pFrac, exp);
        }
      }
    }
  }
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

  int nTime,nPlev,nLat,nLon;
  double radius = 6371000.0;
  double invDlambda = 1.0/(2.0*dlambda);
  double invDphi = 1.0/(2.0*dphi);
  double coef;

  nTime = U->get_dim(0)->size();
  nPlev = U->get_dim(1)->size();
  nLat = U->get_dim(2)->size();
  nLon = U->get_dim(3)->size();


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
          -3.0*UMat[t][p][0][b]*cosphi[0])*invDphi;
        dUdphi[t][p][nLat-1][b]=(3.0*UMat[t][p][nLat-1][b]*cosphi[nLat-1]\
          -4.0*UMat[t][p][nLat-2][b]*cosphi[nLat-2]\
          +UMat[t][p][nLat-3][b]*cosphi[nLat-3])*invDphi;
        for (int a=1; a<(nLat-1); a++){
          dUdphi[t][p][a][b]=(UMat[t][p][a+1][b]*cosphi[a+1]\
           -UMat[t][p][a-1][b]*cosphi[a-1])*invDphi;
        }
      }
    }
  }

  //V WRT LAMBDA
  DataMatrix4D<double> dVdl(nTime, nPlev, nLat, nLon);
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        dVdl[t][p][a][0]=(VMat[t][p][a][1]-VMat[t][p][a][nLon-1])*invDlambda;
        dVdl[t][p][a][nLon-1]=(VMat[t][p][a][0]-VMat[t][p][a][nLon-2])*invDlambda;
        for (int b=1; b<(nLon-1); b++){
          dVdl[t][p][a][b]=(VMat[t][p][a][b+1]-VMat[t][p][a][b-1])*invDlambda;
        }
      }
    }
  }


  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        coef = 1.0/(radius*cosphi[a]);
        for (int b=0; b<nLon; b++){
          RVMat[t][p][a][b] = coef*(dVdl[t][p][a][b]-dUdphi[t][p][a][b]);
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


  int nTime,nPlev,nLat,nLon;
  double invdp,invdp1,invdp2;
  double invdphi= 1.0/(2.0*dphi);
  double invdlambda = 1.0/(2.0*dlambda);
  double radius = 6371000.0;
  double coef1,coef2,corvar;

  nTime = U->get_dim(0)->size();
  nPlev = U->get_dim(1)->size();
  nLat = U->get_dim(2)->size();
  nLon = U->get_dim(3)->size();

  //Input matrices
  DataMatrix4D<double> UMat(nTime, nPlev, nLat, nLon);
  DataMatrix4D<double> VMat(nTime, nPlev, nLat, nLon);

  //Load data
  U->set_cur(0,0,0,0); 
  U->get(&(UMat[0][0][0][0]),nTime, nPlev, nLat, nLon);

  V->set_cur(0,0,0,0);
  V->get(&(VMat[0][0][0][0]), nTime, nPlev, nLat, nLon);

  //Pressure axis values
  DataVector<double> pVec(nPlev);
  pVals->set_cur((long) 0);
  pVals->get(&(pVec[0]), nPlev);

  //Matrices for the partials
  //PT, U, V WRT P
  DataMatrix4D<double> dpt_dp(nTime, nPlev, nLat, nLon);
  DataMatrix4D<double> du_dp(nTime, nPlev, nLat, nLon);
  DataMatrix4D<double> dv_dp(nTime, nPlev, nLat, nLon);

  invdp1 = 1.0/(2.0*std::fabs(pVec[1]-pVec[0]));
  invdp2 = 1.0/(2.0*std::fabs(pVec[nPlev-1]-pVec[nPlev-2]));

  for (int t=0; t<nTime; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
      //0 case
        dpt_dp[t][0][a][b] = (-PTMat[t][2][a][b]+4.0*PTMat[t][1][a][b]-3.0*PTMat[t][0][a][b])*invdp1;
        du_dp[t][0][a][b] = (-UMat[t][2][a][b]+4.0*UMat[t][1][a][b]-3.0*UMat[t][0][a][b])*invdp1;
        dv_dp[t][0][a][b] = (-VMat[t][2][a][b]+4.0*VMat[t][1][a][b]-3.0*VMat[t][0][a][b])*invdp1;           //end case
        dpt_dp[t][nPlev-1][a][b] = (3.0*PTMat[t][nPlev-1][a][b]-4.0*PTMat[t][nPlev-2][a][b]\
          +PTMat[t][nPlev-3][a][b])*invdp2;
        du_dp[t][nPlev-1][a][b] = (3.0*UMat[t][nPlev-1][a][b]-4.0*UMat[t][nPlev-2][a][b]\
          +UMat[t][nPlev-3][a][b])*invdp2;
        dv_dp[t][nPlev-1][a][b] = (3.0*VMat[t][nPlev-1][a][b]-4.0*VMat[t][nPlev-2][a][b]\
          +VMat[t][nPlev-3][a][b])*invdp2;
        for (int p=1; p<(nPlev-1); p++){
          invdp = 1.0/(2.0*std::fabs(pVec[p+1]-pVec[p]));
          dpt_dp[t][p][a][b] = (PTMat[t][p+1][a][b]-PTMat[t][p-1][a][b])*invdp;
          du_dp[t][p][a][b] = (UMat[t][p+1][a][b]-UMat[t][p-1][a][b])*invdp;
          dv_dp[t][p][a][b] = (VMat[t][p+1][a][b]-VMat[t][p-1][a][b])*invdp;
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
          -3.0*PTMat[t][p][0][b])*invdphi;
        dpt_dphi[t][p][nLat-1][b]=(3.0*PTMat[t][p][nLat-1][b]-4.0*PTMat[t][p][nLat-2][b]\
          +PTMat[t][p][nLat-2][b])*invdphi;
        for (int a=1; a<(nLat-1); a++){
          dpt_dphi[t][p][a][b]=(PTMat[t][p][a+1][b]-PTMat[t][p][a-1][b])*invdphi;
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
        dpt_dl[t][p][a][0]=(PTMat[t][p][a][1]-PTMat[t][p][a][nLon-1])*invdlambda;
        dpt_dl[t][p][a][nLon-1]=(PTMat[t][p][a][nLon-2]-PTMat[t][p][a][0])*invdlambda;
        for (int b=1; b<(nLon-1); b++){
          dpt_dl[t][p][a][b]=(PTMat[t][p][a][b+1]-PTMat[t][p][a][b-1])*invdlambda;
        }
      }
    }
  }
  coef2 = 1.0/radius;
  DataMatrix4D<double> PVMat(nTime, nPlev, nLat, nLon);
  //PV Calculation!
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        coef1 = 1.0/(radius*cosphi[a]);
        corvar=coriolis[a];
        for (int b=0; b<nLon; b++){
          PVMat[t][p][a][b] = 9.80616*(coef1*dv_dp[t][p][a][b]*dpt_dl[t][p][a][b]\
           -coef2*du_dp[t][p][a][b]*dpt_dphi[t][p][a][b]\
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


  DataMatrix3D<double> iPVMat(nTime, nLat, nLon);

  double PVbot;
  double PVmid;
  double PVtop;
  
//Eliminate polar regions from calculations
  int i10 = std::fabs(10/lat_res);
  int i171 = std::fabs(171/lat_res);
  double modLevLen = pos_bot-pos_top;
  double invLevLen = 1.0/(2.0*modLevLen);
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
        iPVMat[t][a][b] = (PVtop+PVbot+PVmid)*invLevLen;
      }
    }
  }
 
  intPV->set_cur(0,0,0,0);
  intPV->put(&(iPVMat[0][0][0]),nTime,nLat,nLon);
  std::cout<<"Finished integrating PV."<<std::endl;
} 
//////////////////////////////////////////////////////////////////////////////////////////////////

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
              double PVAnom){

  int nTime,nLat,nLon,nSteps,avgDay,nOutTime;
  double tRes;

  nTime = inIPV->get_dim(0)->size();
  nLat = inIPV->get_dim(1)->size();
  nLon = inIPV->get_dim(2)->size();

//keep sign consistent!
  if (PVAnom<0){
    PVAnom = -PVAnom;
  }

//input IPV
  DataMatrix3D<double> IPVMat(nTime,nLat,nLon);
  inIPV->set_cur(0,0,0);
  inIPV->get(&(IPVMat[0][0][0]),nTime,nLat,nLon);

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

//avg IPV
  avgDay = avgIPV->get_dim(0)->size();
  DataMatrix3D<double> avgMat(avgDay,nLat,nLon);
  avgIPV->set_cur(0,0,0);
  avgIPV->get(&(avgMat[0][0][0]),avgDay,nLat,nLon);

  DataVector<int> avgDayVec(avgDay);
  avgTime->set_cur((long) 0);
  avgTime->get(&(avgDayVec[0]),avgDay);


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
  std::cout<<"There are "<<nDays<<" days in file."<<std::endl;


//Deal with skipped days          
  int d=0;
  DataVector<double> newTime(nOutTime);

  int leapYear=0;
  int leapMonth=0;
  int leapDay=0;
  int leapHour=0;

  for (int t=0; t<nTime; t++){
    if (leap){
      ParseTimeDouble(strTimeUnits, strCalendar, timeVec[t], leapYear,\
        leapMonth, leapDay, leapHour);
      if (leapMonth==2 && leapDay == 29){
        std::cout<<"Leap day! Skipping day."<<std::endl;
        t+=nSteps;
      }
    }
    if (t>=nTime){
      break;
    }
    int nDayIncrease = d/nSteps;
    int currAvgIndex = startAvgIndex + nDayIncrease;
    if (currAvgIndex>364){
      currAvgIndex-=365;
    }
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        devMat[t][a][b] = IPVMat[t][a][b]-avgMat[currAvgIndex][a][b];
      }
    }
    newTime[d] = timeVec[t];
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
  outDev->set_cur(0,0,0);
  outDev->put(&(devMat[0][0][0]),nOutTime,nLat,nLon);
  std::cout<<"Wrote devs to file."<<std::endl;

  outADev->set_cur(0,0,0);
  outADev->put(&(aDevMat[0][0][0]),nOutTime,nLat,nLon);
  std::cout<<"Wrote smoothed devs to file."<<std::endl;
//Divide matrix by PV anomaly value 
//We are looking for negative anomalies in NH and positive anomalies in SH
  DataVector<double> latVec(nLat);
  lat->set_cur((long) 0);
  lat->get(&(latVec[0]),nLat);

  DataMatrix3D<int> posIntDevs(nOutTime,nLat,nLon);
  double invAnom = 1.0/PVAnom;
  double divDev,pos,neg;
  for (int t=0; t<nOutTime; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        divDev = aDevMat[t][a][b]*invAnom;
        //SH: positive anomalies
        if (latVec[a]<0){
          pos = (divDev+std::fabs(divDev))*0.5;
          posIntDevs[t][a][b] = int(pos);
        }
        //NH: negative anomalies
        else if (latVec[a]>=0){
          neg = (divDev-std::fabs(divDev))*0.5;
          posIntDevs[t][a][b] = -int(neg);
        }
      }
    }
  }
  outPosIntDev->set_cur(0,0,0);
  outPosIntDev->put(&(posIntDevs[0][0][0]),nOutTime,nLat,nLon);
  std::cout<<"Wrote integer values to file."<<std::endl;
}


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
      if (a==32 && b== 50){
        std::cout<<"Standard deviation value: "<<std::sqrt(variance)<<std::endl;
      }
    }
  }
}


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
              NcVar *stdDevVar){

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
  std::cout<<"There are "<<nDays<<" days in file."<<std::endl;


//Deal with skipped days          
  int d=0;
  DataVector<double> newTime(nOutTime);

  int leapYear=0;
  int leapMonth=0;
  int leapDay=0;
  int leapHour=0;

  for (int t=0; t<nTime; t++){
    if (leap){
      ParseTimeDouble(strTimeUnits, strCalendar, timeVec[t], leapYear,\
        leapMonth, leapDay, leapHour);
      if (leapMonth==2 && leapDay == 29){
        std::cout<<"Leap day! Skipping day."<<std::endl;
        t+=nSteps;
      }
    }
    if (t>=nTime){
      break;
    }
    int nDayIncrease = d/nSteps;
    int currAvgIndex = startAvgIndex + nDayIncrease;
    if (currAvgIndex>364){
      currAvgIndex-=365;
    }
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        devMat[t][a][b] = GHMat[t][a][b]-avgMat[currAvgIndex][a][b];
      }
    }
    newTime[d] = timeVec[t];
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
        denom = std::fabs(std::sin(latVec[a]*pi/180));
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
//We are looking for negative anomalies in NH and positive anomalies in SH


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
       invAnom = 1./180.;
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
  std::cout<<"ContCheck is " << contCheck << " and tRes is "<<\
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
  std::cout<<"The new date index is "<<dateIndex<<" at the end of the missing fill."<<std::endl;
}

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
