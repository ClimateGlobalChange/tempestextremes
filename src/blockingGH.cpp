///////////////////////////////////////////////////
///          \file blockingGH.cpp
///          \author Marielle Pinheiro
///          \version November 5, 2015

/*This code is an implementation of Tibaldi and Molteni 1990
(generally referred to as TM). A 500-mb geopotential height gradient
is calculated over a specified range of latitudes for a single longitude.
If the two criteria are met for one of 3 dlats, then the longitude is
instantaneously blocked. 
Spatial criterion for block: 12 degrees contiguous latitude blocked 
Temporal criterion for block: 5 days persistent blocking for single latitude

*/

#include "CommandLine.h"
#include "Announce.h"
#include "Exception.h"
#include "NetCDFUtilities.h"
#include "netcdfcpp.h"
#include "DataVector.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"
#include "interp_z500.h"
#include "blockingUtilities.h"

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
    //name of 2d file
    std::string file2d;
    //input variable
    std::string varName;
    //interpolate, yes or no?
    bool interp_check;
    //If Pressure is hPa, check for 500 instead of 50000
    bool is_hPa;
    //is variable 4d?
    bool v4d;


    //Parse command line
    BeginCommandLine()
      CommandLineString(fileIn, "in","")
      CommandLineString(outFile,"out","")
      CommandLineString(file2d,"in2d","")
      CommandLineString(varName,"varname","")
      CommandLineBool(interp_check,"ipl")
      CommandLineBool(is_hPa,"hpa")
      CommandLineBool(v4d,"is4d")
      ParseCommandLine(argc,argv);
    EndCommandLine(argv);
    AnnounceBanner();

    if (fileIn == ""){
      _EXCEPTIONT("No input file (--in) specified");
    }
    if (outFile == ""){
      std::string fileInCopy = fileIn;
      outFile = fileInCopy.replace(fileInCopy.end()-3,fileInCopy.end(),"_GHblock.nc");
    }

    //If variable needs interpolation, do it first
    if (interp_check){
      if (file2d ==""){
        _EXCEPTIONT("No input file (--in2d) specified");
      }
      NcFile readin_int(fileIn.c_str());
      if (!readin_int.is_valid()){
        _EXCEPTION1("Unable to open file \"%s\" for reading",
          fileIn.c_str());
      }
      //z500 file name
      std::string zfile = fileIn.replace(fileIn.end()-3,fileIn.end(),"_z500.nc");
      NcFile readin_out(zfile.c_str(),NcFile::Replace, NULL,0,NcFile::Offset64Bits);
      if (!readin_out.is_valid()){
        _EXCEPTION1("Unable to open file \"%s\" for writing",zfile.c_str());
      }
      //interpolate file
      interp_z500(readin_int, file2d, readin_out);
      readin_out.close();
      fileIn = zfile;

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

    if (varName == "Z"){
  //Variables for calculations
      NcVar *zvar = readin.get_var("Z");

  //If the variable is 4D, need to get only 500 mb data
      if (v4d){
        NcDim *lev = readin.get_dim("lev");
        int nLev = lev->size();
        NcVar *levvar = readin.get_var("lev");

        DataVector<double> pVec(nLev);
        levvar->set_cur((long) 0);
        levvar->get(&(pVec[0]),nLev);

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
      }
      //otherwise, load the 3D data
      else{
        zvar->set_cur(0,0,0);
        zvar->get(&(ZData[0][0][0]),nTime,nLat,nLon);
      }

      //NOTE: Geopotential needs to be divided by g in order to get height

      for (int t=0;t<nTime;t++){
        for (int a=0; a<nLat;a++){
          for (int b=0; b<nLon; b++){
            ZData[t][a][b]/=9.8;
          }
        }
      }
    }
    else if (varName == "GH"){
      NcVar *zvar = readin.get_var("GH");
      zvar->set_cur(0,0,0);
      zvar->get(&(ZData[0][0][0]),nTime,nLat,nLon);

    }
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
    if (time_vals->get_att("calendar") == NULL){
      time_vals->add_att("calendar","standard");
    }
    copy_dim_var(latvar, lat_vals);
    copy_dim_var(lonvar, lon_vals);

  //FOR DEBUG: WRITE 500 mb GH to file
  NcVar *gh_vals = file_out.add_var("GH",ncDouble,out_time,out_lat,out_lon);
  gh_vals->set_cur(0,0,0);
  gh_vals->put(&(ZData[0][0][0]),nTime,nLat,nLon); 

  //get the spatial resolution 
    DataVector<double> latVec(nLat);
    latvar->set_cur((long) 0);
    latvar->get(&(latVec[0]),nLat);


    double latRes = latVec[0]-latVec[1];
    std::cout<<"latres is "<<latRes<<std::endl;
    double latDelta = 15./latRes;
    int iDelta = (int) latDelta;
    std::cout<<"delta is "<<iDelta<<std::endl;

  //Calculate blocking index for 75-15N, 15-75S
    int NHLatStart;
    int NHLatEnd;
    int SHLatStart;
    int SHLatEnd;

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
    std::cout<<"NH Start/End: "<<NHLatStart<<"/"<<NHLatEnd<<std::endl;
    std::cout<<"SH Start/End: "<<SHLatStart<<"/"<<SHLatEnd<<std::endl;
  
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
