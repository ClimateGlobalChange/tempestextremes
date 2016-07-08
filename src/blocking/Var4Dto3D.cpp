////////////////////////////////////////////////
///         \file var4Dto3D.cpp
///         \author Marielle Pinheiro
///         \version December 2, 2015

/*This code takes a 4D variable and outputs a 3D variable (at the moment
only along the lev axis)*/

#include "CommandLine.h"
#include "Announce.h"
#include "Exception.h"
#include "NetCDFUtilities.h"
#include "netcdfcpp.h"
#include "DataVector.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"
#include "BlockingUtilities.h"

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>

int main(int argc, char** argv){
  NcError error(NcError::verbose_nonfatal);
  try{

    std::string fileIn;
    std::string fileOut;
    bool is_hPa;

    BeginCommandLine()
      CommandLineString(fileIn,"in","");
      CommandLineString(fileOut,"out","");
      CommandLineBool(is_hPa,"hpa");
      ParseCommandLine(argc,argv);
    EndCommandLine(argv);

    if (fileIn == ""){
      _EXCEPTIONT("No input file (--in) specified");
    }

    if (fileOut == ""){
      std::string fileInCopy = fileIn;
      fileOut = fileInCopy.replace(fileInCopy.end()-3,fileInCopy.end(),"_GH.nc");
    }

    NcFile readin(fileIn.c_str());

    //Dimensions and associated variables
    NcDim *time = readin.get_dim("time");
    int nTime = time->size();
    NcVar *timevar = readin.get_var("time");

    NcDim *lev = readin.get_dim("lev");
    int nLev = lev->size();
    NcVar *levvar = readin.get_var("lev");

    DataVector<double> pVec(nLev);
    levvar->set_cur((long) 0);
    levvar->get(&(pVec[0]),nLev);

    NcDim *lat = readin.get_dim("lat");
    int nLat = lat->size();
    NcVar *latvar = readin.get_var("lat");

    NcDim *lon = readin.get_dim("lon");
    int nLon = lon->size();
    NcVar *lonvar = readin.get_var("lon");

    //Variables for calculations
    NcVar *zvar = readin.get_var("Z");
    DataMatrix3D<double> ZData(nTime, nLat, nLon);

    double pval = 50000.0;
    if (is_hPa){
      pval = 500.0;
    }
    int pIndex;
    for (int x=0; x<nLev; x++){
      if (std::fabs(pVec[x]-pval)<0.0001){
        pIndex = x;
        break;
      }
    }
    std::cout<<"pIndex: "<<pIndex<<std::endl;
    zvar->set_cur(0,pIndex,0,0);
    zvar->get(&(ZData[0][0][0]),nTime,1,nLat,nLon);

    for (int t=0; t<nTime; t++){
      for (int a=0; a<nLat;a++){
        for (int b=0; b<nLon; b++){
          ZData[t][a][b]/=9.8;
        }
      }
    }
    NcFile file_out(fileOut.c_str(),NcFile::Replace,NULL,0,NcFile::Offset64Bits);

    //Dimensions: time, lat, lon
    NcDim *out_time = file_out.add_dim("time", nTime);
    NcDim *out_lat = file_out.add_dim("lat", nLat);
    NcDim *out_lon = file_out.add_dim("lon", nLon);

  //COPY EXISTING DIMENSION VALUES
    NcVar *time_vals = file_out.add_var("time", ncDouble, out_time);
    NcVar *lat_vals = file_out.add_var("lat", ncDouble, out_lat);
    NcVar *lon_vals = file_out.add_var("lon", ncDouble, out_lon);

    copy_dim_var(timevar, time_vals);
    if (time_vals->get_att("calendar") == NULL){
      time_vals->add_att("calendar","standard");
    }
    copy_dim_var(latvar, lat_vals);
    copy_dim_var(lonvar, lon_vals);

    NcVar *gh_vals = file_out.add_var("GH",ncDouble,out_time,out_lat,out_lon);
    gh_vals->set_cur(0,0,0);
    gh_vals->put(&(ZData[0][0][0]),nTime,nLat,nLon);

    file_out.close();
  }catch (Exception & e){
    std::cout<<e.ToString() <<std::endl;
  }
}
