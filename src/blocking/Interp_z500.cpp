//////////////////////////////////
///
///    \file interpolate.cpp
///    \author Marielle Pinheiro
///    \version March 24, 2015

#include "Interp_z500.h"
#include "BlockingUtilities.h"
#include "NetCDFUtilities.h"
#include "netcdfcpp.h"
#include "DataVector.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"
#include <cstdlib>
#include <cmath>
#include <cstring>

//////////////////////////////////////////////////////////////////////// 

void interp_1lev(NcVar *var,
                     NcVar *hyam,
                     NcVar *hybm,
                     NcVar *ps,
                     double plev,
                     NcVar *NewVar
){

  int nTime,nLev,nLat,nLon;
  double A1,A2,B1,B2,p1,p2,weight;

  //Array dimensions
  nTime = var->get_dim(0)->size();
  nLev = var->get_dim(1)->size();
  nLat = var->get_dim(2)->size();
  nLon = var->get_dim(3)->size();

 //Matrix to store PS
  DataMatrix3D <double> matPS(nTime, nLat, nLon);
  ps->set_cur(0, 0, 0);
  ps->get(&(matPS[0][0][0]), nTime, nLat, nLon);

  //Matrix to store input variable data
  DataMatrix4D<double> matVar(nTime, nLev, nLat, nLon);
  var->set_cur(0, 0, 0, 0);
  var->get(&(matVar[0][0][0][0]), nTime, nLev, nLat, nLon);

  //Matrix to store output variable data
  DataMatrix3D<double> matVarOut(nTime, nLat, nLon);

  //hybrid coefficient A
  DataVector<double> vecHyam(nLev);
  hyam->set_cur((long) 0);
  hyam->get(&(vecHyam[0]), nLev);

  //hybrid coefficient B
  DataVector<double> vecHybm(nLev);
  hybm->set_cur((long) 0);
  hybm->get(&(vecHybm[0]), nLev);

  //Loop over input data and interpolate to output var
  for (int t=0; t<nTime; t++){
    for (int l=0; l<(nLev-1); l++){
      A1 = vecHyam[l];
      B1 = vecHybm[l];
      A2 = vecHyam[l+1];
      B2 = vecHybm[l+1];
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          p1 = 100000.0 * A1 + matPS[t][a][b] * B1;
          p2 = 100000.0 * A2 + matPS[t][a][b] * B2;
          if (p1<plev && p2>=plev){
              weight = ((plev-p1)/(p2-p1));
              matVarOut[t][a][b] = weight*matVar[t][l+1][a][b]
                               + (1.0-weight)*matVar[t][l][a][b];
          }
        }
      }
    }
  }
  std::cout<<"Finished interpolating variable.\n";

  NewVar->set_cur(0, 0, 0);
  NewVar->put(&(matVarOut[0][0][0]), nTime, nLat, nLon);
  CopyNcVarAttributes(var, NewVar);
}



void interp_z500(NcFile & readin, 
                 const std::string & strname_2d, 
                 const std::string & varname,
                 NcFile & ifile_out) {

  //open 2D PS file
  NcFile readin_2d(strname_2d.c_str());
  if (!readin_2d.is_valid()) {
    _EXCEPTION1("Unable to open file \"%s\" for reading",
      strname_2d.c_str());
  }

  //Dimensions and associated variables
  NcDim *time = readin.get_dim("time");
  int time_len = time->size();
  NcVar *timevar = readin.get_var("time");

  NcDim *lev = readin.get_dim("lev");
  int lev_len = lev->size();

  NcDim *lat = readin.get_dim("lat");
  int lat_len = lat->size();
  NcVar *latvar = readin.get_var("lat");

  NcDim *lon = readin.get_dim("lon");
  int lon_len = lon->size();
  NcVar *lonvar = readin.get_var("lon");

  //Variables
  NcVar *zvar = readin.get_var(varname.c_str());
  
  //2D variables
  NcVar *ps = readin_2d.get_var("PS");
  NcVar *hyam = readin.get_var("hyam");
  NcVar *hybm = readin.get_var("hybm");
    
  //Write information to outfile
  NcDim *itime = ifile_out.add_dim("time", time_len);
  NcDim *ilat = ifile_out.add_dim("lat", lat_len);
  NcDim *ilon = ifile_out.add_dim("lon", lon_len);
    
  NcVar *itime_vals = ifile_out.add_var("time", ncDouble, itime);
  NcVar *ilat_vals = ifile_out.add_var("lat", ncDouble, ilat);
  NcVar *ilon_vals = ifile_out.add_var("lon", ncDouble, ilon);

  copy_dim_var(timevar, itime_vals);
  copy_dim_var(latvar, ilat_vals);
  copy_dim_var(lonvar, ilon_vals);

  //Add interpolated variables to interpolated outfile
  NcVar *iz = ifile_out.add_var("Z500",, ncDouble, itime, ilat, ilon);
  interp_1lev(zvar, hyam, hybm, ps, 50000.0, iz);


  readin_2d.close();  
  std::cout<<"Finished interpolating Z to 500 mb level."<<std::endl;
} 
