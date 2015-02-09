/////////////////////////////////////////
///
///    \file CLIVAR_blocks.cpp
//     \author Marielle Pinheiro
///    \version February 3, 2015
///


#include "CLIVAR_blocks.h"
#include "Exception.h"
#include "NetCDFUtilities.h"
#include "netcdfcpp.h"

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>


NcVar *interpolate_lev(NcVar *var, NcVar *hyam, NcVar *hybm, NcVar *ps, NcFile outfile){
  //need to calculate at each level, interpolate to new variable
  //How to get dimensions? Ask!

  //Create new variable and copy time, lat, lon dimensions
  //but initialize pressure levels
  //Possible to copy dimensions?
  //Or should the variable be stored as an array, 
  //then written as a new variable?
  NcAtt *att_units = var->get_att("units");
  char *units = att_units->as_string(0);
  NcVar *NewVar = outfile.add_var(units.c_str());

  int nTime = var->get_dim(0)->size();
  int nLev = var->get_dim(1)->size();
  int nLat = var->get_dim(2)->size();
  int nLon = var->get_dim(3)->size();

  DataMatrix3D<double> matPS(nTime, nLat, nLon);
  var->set_cur(0, 0, 0);
  var->get(&(matPS[0][0][0]), nTime, nLat, nLon);

  DataMatrix4D<double> matVar(nTime, nLev, nLat, nLon);
  var->set_cur(0, 0, 0, 0);
  var->get(&(matVar[0][0][0][0]), nTime, nLev, nLat, nLon);

  DataMatrix4D<double> matVarOut(nTime, nLev, nLat, nLon);

  DataVector<float> vecHyam(nLev);
  hyam->set_cur((long) 0);
  hyam->get(&(vecHyam[0]), nLev);

  DataVector<float> vecHybm(nLev);
  hybm->set_cur((long) 0);
  hybm->get(&(vecHybm[0]), nLev);
  
  //Note: need to define these lengths!
  for (t=0; t++; t<time_len){
    for (l=0; l++; l<(lev_len-1)){
      float A1 = vecHyam[l];
      float B1 = vecHybm[l];
      float A2 = vecHyam[l+1];
      float B2 = vecHybm[l+1];
      for (a=0; a++; a<lat_len){
        for (b=0; b++; b<lon_len){
          float p1 = 100000.0 * A1 + matPS[t][a][b] * B1;
          float p2 = 100000.0 * A2 + matPS[t][a][b] * B2;
          for (p=0; p++; p<plev_len){
            if (p1<plev[p] && p2>=plev[p]){
              float weight = ((plev[p]-p1)/(p2-p1));
              matVarOut[t][p][a][b] = weight*matVar[t][l+1][a][b]
                               + (1.-weight)*matVar[t][l][a][b];
            }
          }
        }
      }
    }
  }

  NewVar->set_cur(0, 0, 0, 0);
  NewVar->put(&(matVarOut[0][0][0][0]), nTime, nLev, nLat, nLon);

  return(NewVar);
}
//////////////////////////////////////////////////////////////////////// 
int main(int argc, char **argv){
  if (argc<3){
    cout<< "Error: provide input file names for both 3D and 2D.";
    return 1;
  }
  std::string strfile_in = argv[1];

  //Open input file (as read-only)
  NcFile readin(strfile_in.c_str());
  cout<< "Reading file "<< strfile_in.c_str();  

  //Open 2D input file
  std::string strfile_2d = argv[2];
  NcFile readin_2d(strfile_2d.c_str());
  cout<< "Reading file " << strfile_2d.c_str();
  
  //Dimensions
  //What does the arrow symbol do?
  NcDim *time = readin.get_dim("time");
  long time_len = time->size();

  NcDim *lev = readin.get_dim("lev");
  long lev_len = lev->size();

  NcDim *lat = readin.get_dim("lat");
  long lat_len = lat->size();

  NcDim *lon = readin.get_dim("lon");
  long lon_len = lon->size();

  //Coefficients
  NcVar *hyam = readin.get_var("hyam");
  NcVar *hybm = readin.get_var("hybm");

  //Variables
  NcVar *temp = readin.get_var("T");
  NcVar *uvar = readin.get_var("U");
  NcVar *vvar = readin.get_var("V");
  NcVar *ps = readin_2d.get_var("PS");

  //Create output file
  std::string strfile_out = strfile_in.replace(".nc", "test.nc");
  NcFile file_out(strfile_out.c_str(), NcFile::Replace);
  NcDim *out_time = file_out.add_dim("time", time_len);
  //How to create pressure level?
  NcDim *out_plev = file_out.add_dim("lev", plev_len);
  NcDim *out_lat = file_out.add_dim("lat", lat_len);
  NcDim *out_lon = file_out.add_dim("lon", lon_len);
  

  ////////////////////////////////////////////////////////////
  //Close input files
  filein.close()
  filein_2d.close()
}
