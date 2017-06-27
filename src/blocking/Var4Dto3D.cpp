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
#include "Interpolate.h"

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>

int main(int argc, char** argv){
  NcError error(NcError::verbose_nonfatal);
  try{

    std::string fileIn;
    std::string fileIn_2D;
    std::string fileOut;
    bool is_hPa;
    bool interp_check;
    bool ZtoGH;
    std::string varlist;

    BeginCommandLine()
      CommandLineString(fileIn,"in","");
      CommandLineString(fileIn_2D,"in2d","");
      CommandLineString(fileOut,"out","");
      CommandLineBool(is_hPa,"hpa");
      CommandLineBool(interp_check,"ipl");
      CommandLineBool(ZtoGH,"gh");
      CommandLineString(varlist,"varlist","");
      ParseCommandLine(argc,argv);
    EndCommandLine(argv);

    if (fileIn == ""){
      _EXCEPTIONT("No input file (--in) specified");
    }

    if (fileOut == ""){
      std::string fileInCopy = fileIn;
      fileOut = fileInCopy.replace(fileInCopy.end()-3,fileInCopy.end(),"_z500.nc");
    }
    if (varlist == ""){
       _EXCEPTIONT("Need to provide variable names with --varlist flag.");
    }
    //if variable needs to be interpolated, do that first!
    if (interp_check){
      if (fileIn_2D==""){
        _EXCEPTIONT("No input file (--in2d) specified for surface variables");
      }
      if (varlist==""){
        _EXCEPTIONT("Need to provide at least 1 variable name with the --varlist flag");
      }
      NcFile interp_in(fileIn.c_str());
      std::string interp_outname = fileIn.replace(fileIn.end()-3,fileIn.end(),"_ipl_3D.nc");
      //open output file for interpolated variable
      NcFile interp_out(interp_outname.c_str(),NcFile::Replace, NULL, 0, NcFile::Offset64Bits);
      if (!interp_out.is_valid()){
        _EXCEPTIONT("Unable to open file for interpolated variables");
      }
      //Interpolate variables specified by list into file
      interp_util(interp_in,fileIn_2D,varlist,interp_out);
      interp_out.close();
      //set input file name to the name of the file that was interpolated
      fileIn = interp_outname;

    }

  
    NcFile readin(fileIn.c_str());

    //Dimensions and associated variables
    NcDim *time = readin.get_dim("time");
    int nTime = time->size();
    NcVar *timevar = readin.get_var("time");

    NcDim *lev = readin.get_dim("lev");
    int nLev = lev->size();
    NcVar *levvar = readin.get_var("lev");

    //Create a data vector with the associated pressure values 
    DataVector<double> pVec(nLev);
    levvar->set_cur((long) 0);
    levvar->get(&(pVec[0]),nLev);

    NcDim *lat = readin.get_dim("lat");
    int nLat = lat->size();
    NcVar *latvar = readin.get_var("lat");

    NcDim *lon = readin.get_dim("lon");
    int nLon = lon->size();
    NcVar *lonvar = readin.get_var("lon");

    //Find the index of the 500 mb level
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

    //Open output file
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

    //Split var list  
    std::string delim = ",";
    size_t pos = 0;
    std::string token;
    std::vector<std::string> varVec;

    while((pos = varlist.find(delim)) != std::string::npos){
      token = varlist.substr(0,pos);
      varVec.push_back(token);
      varlist.erase(0,pos + delim.length());
    }
    varVec.push_back(varlist);
    for (int v=0; v<varVec.size(); v++){
      std::cout<<"Vector contains string "<<varVec[v].c_str()<<std::endl;

      NcVar *vvar = readin.get_var(varVec[v].c_str());
      DataMatrix3D<double> VData(nTime, nLat, nLon);

      vvar->set_cur(0,pIndex,0,0);
      vvar->get(&(VData[0][0][0]),nTime,1,nLat,nLon);
      if (varVec[v]=="Z" && ZtoGH){
        for (int t=0; t<nTime; t++){
          for (int a=0; a<nLat; a++){
            for (int b=0; b<nLon; b++){
              VData[t][a][b]/=9.8;
            }
          }
        }
      }

      NcVar *outvar = file_out.add_var(varVec[v].c_str(),ncDouble,out_time,out_lat,out_lon);
      outvar->set_cur(0,0,0);
      outvar->put(&(VData[0][0][0]),nTime,nLat,nLon);
    }


/*
    for (int t=0; t<nTime; t++){
      for (int a=0; a<nLat;a++){
        for (int b=0; b<nLon; b++){
          ZData[t][a][b]/=9.8;
        }
      }
    }
*/


    file_out.close();
  }catch (Exception & e){
    std::cout<<e.ToString() <<std::endl;
  }
}
