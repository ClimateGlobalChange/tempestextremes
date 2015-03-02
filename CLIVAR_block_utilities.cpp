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
){
  //Array dimensions
  int nTime = var->get_dim(0)->size();
  int nLev = var->get_dim(1)->size();
  int nLat = var->get_dim(2)->size();
  int nLon = var->get_dim(3)->size();
  int npLev = pLev->get_dim(0)->size();

  std::cout<<"Dimensions: Time: "<<nTime<<", hybrid level: "<<nLev<< " lat: " << nLat << " lon: " << nLon << " plev: " << npLev<<std::endl;

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
//          if (interp_check == 0){
//            std::cout<< "Did not interpolate this variable to pressure levels for (l,a,b) = ("<< l <<", "<<a<<", "<<b<<")"<<std::endl;
//          }
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

////////////////////////////////////////////

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
){ 
  double radius = 6371000.0;
  double pi = 4.0*std::atan(1.0);
  double radian = 180.0/pi;
  double sigma = std::pow(7.2921, -5.0);

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

  lat_res = latVec[1]-latVec[0];
  lon_res = lonVec[1]-lonVec[0];
  p_res = pVec[1]-pVec[0];
 
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
	NcVar *PT
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
  DataMatrix4D<double> PTMat(nTime, nPlev, nLat, nLon);
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
  PT->set_cur(0,0,0,0);
  PT->put(&(PTMat[0][0][0][0]), nTime, nPlev, nLat, nLon);
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
	NcVar *vorticity
){
  std::cout<<"Starting rel vort calculations."<<std::endl;

  int nTime = U->get_dim(0)->size();
  int nPlev = U->get_dim(1)->size();
  int nLat = U->get_dim(2)->size();
  int nLon = U->get_dim(3)->size();

 std::cout<<"Dimensions: Time: "<<nTime<<", lat: " << nLat << " lon: " << nLon << " plev: " << nPlev<<std::endl;
 double radius = 6371000.0;
  std::cout<<"Initializing data matrices."<<std::endl;
//INPUTS
  DataMatrix4D<double> VMat(nTime, nPlev, nLat, nLon);

  std::cout<< "About to load variables."<<std::endl;
  V->set_cur(0,0,0,0);
  std::cout<<"Set V to 0."<<std::endl;
  V->get(&(VMat[0][0][0][0]), nTime, nPlev, nLat, nLon);
  std::cout<<"Loaded V matrix."<<std::endl;

  DataMatrix4D<double> UMat(nTime, nPlev, nLat, nLon);
  U->set_cur(0,0,0,0);
  U->get(&(UMat[0][0][0][0]),nTime, nPlev, nLat, nLon);
  std::cout<<"Set U data matrix."<<std::endl;
  std::cout<<"Initializing output data matrix."<<std::endl;
//OUTPUT: Relative vorticity
  DataMatrix4D<double> RVMat(nTime, nPlev, nLat, nLon);
  for (int t=0; t<nTime; t++){
    for (int p=0; p<nPlev; p++){
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          double dvdlambda;
          double dudphi;
	  //END CASES (PERIODIC BOUNDARY CONDITION)
          if (b==0){
	    dvdlambda = (VMat[t][p][a][1]-VMat[t][p][a][nLon-1])/(2.0*dlambda);
          }
	  else if (b==(nLon-1)){
	    dvdlambda = (VMat[t][p][a][0]-VMat[t][p][a][nLon-2])/(2.0*dlambda);
	  }
          //OTHERWISE:
          else{
	    dvdlambda = (VMat[t][p][a][b+1]-VMat[t][p][a][b-1])/(2.0*dlambda);
	  }
          //END CASES (POLES)
          if (a==0){
	    dudphi = (-UMat[t][p][2][b]*cosphi[2] + 4.0*UMat[t][p][1][b]*cosphi[1]-3.0*UMat[t][p][0][b]*cosphi[0])/(2.0*dphi);
	  }
	  else if (a==(nLat-1)){
	    dudphi = (3.0*UMat[t][p][nLat-1][b]*cosphi[nLat-1] - 4.0*UMat[t][p][nLat-2][b]*cosphi[nLat-2] + UMat[t][p][nLat-3][b]*cosphi[nLat-3])/(2.0*dphi);
	  }
	  //OTHERWISE:
	  else{
	    dudphi = (UMat[t][p][a+1][b]*cosphi[a+1]-UMat[t][p][a-1][b]*cosphi[a-1])/(2.0*dphi);
	  }
	  RVMat[t][p][a][b] = (1.0/(radius*cosphi[a]))*(dvdlambda-dudphi);
          if (t==0){std::cout<<"Rel vort is "<< RVMat[t][p][a][b]<<std::endl;}
        }
      }
    }
  }

  std::cout<<"Finished calculating relative vorticity."<<std::endl;
  vorticity->set_cur(0,0,0,0);
  vorticity->put(&(RVMat[0][0][0][0]),nTime, nPlev, nLat, nLon);
}

///////////////////////////////////////////////////////////////////////////////////////////////

void PV_calc(
	NcVar *U, 
	NcVar *V, 
	NcVar *PT, 
	NcVar *rVort, 	
	DataVector<double> coriolis,
	double dphi,
	double dlambda,
	double dp,
	NcVar *PV){

  int nTime = U->get_dim(0)->size();
  int nPlev = U->get_dim(1)->size();
  int nLat = U->get_dim(2)->size();
  int nLon = U->get_dim(3)->size();

  //Input matrices
  DataMatrix4D<double> UMat;
  DataMatrix4D<double> VMat;
  DataMatrix4D<double> PTMat;
  DataMatrix4D<double> RVMat;

  //Load data
  U->set_cur(0,0,0,0); 
  U->get(&(UMat[0][0][0][0]),nTime, nPlev, nLat, nLon);

  V->set_cur(0,0,0,0);
  V->get(&(VMat[0][0][0][0]), nTime, nPlev, nLat, nLon);

  PT->set_cur(0,0,0,0);
  PT->get(&(PTMat[0][0][0][0]), nTime, nPlev, nLat, nLon);

  rVort->set_cur(0,0,0,0);
  rVort->get(&(RVMat[0][0][0][0]), nTime, nPlev, nLat, nLon);

  


} 



