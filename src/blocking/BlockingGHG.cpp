///////////////////////////////////////////////////
///          \file blockingGHG.cpp
///          \author Marielle Pinheiro
///          \version February 25, 2017

/*This code is an implementation of Tibaldi and Molteni 1990
(generally referred to as TM). A 500-mb geopotential height gradient
is calculated over a specified range of latitudes for a single longitude.
If the two criteria are met for one of 3 dlats, then the longitude is
instantaneously blocked. 
Spatial criterion for block: 12 degrees contiguous latitude blocked 
Temporal criterion for block: 5 days persistent blocking for single latitude

If data is not in the correct file format, first run Var4Dto3D to generate 
the appropriate geopotential height file
*/

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
    //file which contains Z500
    std::string fileIn;
    //name of output file
    std::string outFile;
    //input variable
    std::string varName;


    BeginCommandLine()
      CommandLineString(fileIn, "in","")
      CommandLineString(outFile,"out","")
      CommandLineString(varName,"varname","")
    ParseCommandLine(argc,argv);
    EndCommandLine(argv);
    AnnounceBanner();

    if (fileIn == ""){
      _EXCEPTIONT("No input file (--in) specified");
    }
    if (varName == ""){
      _EXCEPTIONT("No variable name (--varname) specified");
    }

    if (outFile == ""){
      std::string fileInCopy = fileIn;
      outFile = fileInCopy.replace(fileInCopy.end()-3,fileInCopy.end(),"_GHG.nc");
    }
    NcFile readin(fileIn.c_str());
    if (!readin.is_valid()){
      _EXCEPTION1("Invalid file \"%s\"", fileIn.c_str());
    }
    //Dimensions and associated variables
    NcDim *time = readin.get_dim("time");
    int nTime = time->size();
    NcVar *timevar = readin.get_var("time");

    NcDim *lat = readin.get_dim("lat");
    int nLat = lat->size();
    NcVar *latvar = readin.get_var("lat");

    NcDim *lon = readin.get_dim("lon");
    int nLon = lon->size();
    NcVar *lonvar = readin.get_var("lon");

    DataMatrix3D<double>ZData(nTime,nLat,nLon);
    NcVar *zvar = readin.get_var(varName.c_str());
    zvar->set_cur(0,0,0);
    zvar->get(&(ZData[0][0][0]),nTime,nLat,nLon);

    //create output file
    NcFile file_out(outFile.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);

    //Dimensions: time, lat, lon
    NcDim *out_time = file_out.add_dim("time", nTime);
    NcDim *out_lat = file_out.add_dim("lat", nLat);
    NcDim *out_lon = file_out.add_dim("lon", nLon);

    //COPY EXISTING DIMENSION VALUES
    NcVar *time_vals = file_out.add_var("time", ncDouble, out_time);
    NcVar *lat_vals = file_out.add_var("lat", ncDouble, out_lat);
    NcVar *lon_vals = file_out.add_var("lon", ncDouble, out_lon);

    copy_dim_var(timevar, time_vals);
    copy_dim_var(latvar, lat_vals);
    copy_dim_var(lonvar, lon_vals);

    DataVector<double> latVec(nLat);
    latvar->set_cur((long) 0);
    latvar->get(&(latVec[0]),nLat);


    double latRes = latVec[0]-latVec[1];
    std::cout<<"latres is "<<latRes<<std::endl;
    double latDelta = 15./latRes;
    int iDelta = (int) latDelta;
    std::cout<<"delta is "<<iDelta<<std::endl;

    int NHLatStart;
    int NHLatEnd;
    int SHLatStart;
    int SHLatEnd;

    //Checking blocking for 
    for (int x=0; x<nLat; x++){
      if (std::fabs(75.-latVec[x])<0.0001){
        NHLatStart = x;
      }
      if (std::fabs(35.-latVec[x])<0.0001){
        NHLatEnd = x;
      }
      if (std::fabs(-35-latVec[x])<0.0001){
        SHLatStart = x;
      }
      if (std::fabs(-75.-latVec[x])<0.0001){
        SHLatEnd = x;
      }
    }

  //Make sure that the indices are in the correct order!
    if (NHLatStart>NHLatEnd){
      std::cout<<"Switching NH indices!"<<std::endl;
      int Ntemp = NHLatStart;
      NHLatStart = NHLatEnd;
      NHLatEnd = Ntemp;
    }

    if (SHLatStart>SHLatEnd){
      std::cout<<"Switching SH indices!"<<std::endl;
      int Stemp = SHLatStart;
      SHLatStart = SHLatEnd;
      SHLatEnd = Stemp;
   }

    double z_N, z_S, z_C;
    double lat_N, lat_S, lat_C;
    int i_N, i_S;

  //NEW VARIABLE: TM BLOCKING INDEX
    DataMatrix3D<double> outIndex(nTime,nLat,nLon);
    for (int t=0; t<nTime; t++){
      for (int b=0; b<nLon; b++){
        //Calculate blocking index for NH
        for (int a=NHLatStart; a<=NHLatEnd; a++){
          i_N = a-iDelta;
          i_S = a+iDelta;
          //center values
          lat_C = latVec[a];
          z_C = ZData[t][a][b];
          //upper values
          lat_N = latVec[i_N];
          z_N = ZData[t][i_N][b];
          //lower values
          lat_S = latVec[i_S];
          z_S = ZData[t][i_S][b];

          outIndex[t][a][b] = GHcheck(z_C,z_N,z_S,lat_C,lat_N,lat_S,"N");
        }
        for (int a=SHLatStart; a<=SHLatEnd; a++){
          i_N = a-iDelta;
          i_S = a+iDelta;
          //center values
          lat_C = latVec[a];
          z_C = ZData[t][a][b];
          //upper values
          lat_N = latVec[i_N];
          z_N = ZData[t][i_N][b];
          //lower values
          lat_S = latVec[i_S];
          z_S = ZData[t][i_S][b];

          outIndex[t][a][b] = GHcheck(z_C,z_N,z_S,lat_C,lat_N,lat_S,"S");
        }
      }
    }

  //ADD NEW BLOCKING INDEX
    NcVar *blocking_vals = file_out.add_var("GHGrad", ncDouble, out_time, out_lat, out_lon);
    blocking_vals->set_cur(0,0,0);
    blocking_vals->put(&(outIndex[0][0][0]),nTime,nLat,nLon);

    file_out.close();

  }catch (Exception & e){
    std::cout<<e.ToString() <<std::endl;
  }


}
