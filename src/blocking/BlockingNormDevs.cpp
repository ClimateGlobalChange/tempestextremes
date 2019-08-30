/////////////////////////////////////////////////////////
///
///           \file blockingNormDevs.cpp
///           \author Marielle Pinheiro
///           \version June 1, 2015

/*This code is the third step in the potential vorticity code
based on the Schwierz et al 2004 paper. It calculates deviations 
from the average produced in the previous step (DIPV) and filters them with
2-day smoothing (ADIPV), then divides ADIPV by the specified deviation 
and outputs integer values (INT_ADIPV) which can then be used by StitchBlobs

*/

#include "BlockingUtilities.h"
#include "CommandLine.h"
#include "Announce.h"
#include "Exception.h"
#include "NetCDFUtilities.h"
#include "netcdfcpp.h"
#include "DataArray1D.h"
#include "DataArray3D.h"
#include "DataArray4D.h"
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>


int main(int argc, char **argv){
//  NcError error(NcError::verbose_nonfatal);
  NcError error(NcError::silent_nonfatal);
  try{

    //list of input files for which to calculate normalized deviations
    std::string fileList;
    std::string fileName;
    std::string outputName;
    std::string varName;
    std::string threshName;
    std::string threshVarName;
    std::string tname,latname,lonname;
    std::string insuff,outsuff;
    std::string normName;
    bool PVCalc;
    bool GHCalc;
    bool const_thresh;
    bool appendBool;
    double anomVal,minThresh;
    int startday, endday;
    size_t pos,len;

    BeginCommandLine()
      CommandLineString(fileName,"in","");
      CommandLineString(outputName,"out","");
      CommandLineString(fileList, "inlist", "");
      CommandLineString(varName,"varname","");
//      CommandLineString(avgName, "avg", "");
//      CommandLineString(avgVarName, "avgname","");
      CommandLineString(threshName,"thresh","");
      CommandLineString(threshVarName,"threshname","");
      CommandLineString(normName,"normname","");
      CommandLineString(insuff,"insuff",".nc");
      CommandLineString(outsuff,"outsuff","_norm.nc");
      CommandLineBool(PVCalc,"pv");
      CommandLineBool(GHCalc,"z500");
      CommandLineBool(const_thresh,"const");
      CommandLineBool(appendBool,"append");
      CommandLineDouble(anomVal,"threshval",0.);
      CommandLineDouble(minThresh,"minthreshval",0.);
      CommandLineString(tname,"tname","time");
      CommandLineString(latname,"latname","lat");
      CommandLineString(lonname,"lonname","lon");
      CommandLineInt(startday,"startday",1);
      CommandLineInt(endday,"endday",365);
      ParseCommandLine(argc, argv);
    EndCommandLine(argv);
    AnnounceBanner();

  /*  if ((!PVCalc) && (!GHCalc)){
      _EXCEPTIONT("Need to specify either PV (--pv) or GH (--gh) calculations.");
    }*/
    if (fileName == "" && fileList == ""){
      _EXCEPTIONT("Need to specify either input file (--in) or file list (--inlist).");
    }
    if (fileName != "" && fileList != ""){
      _EXCEPTIONT("Cannot specify both file name (--in) and file list (--inlist).");
    }
    if (outputName != "" && fileList != ""){
      _EXCEPTIONT("Currently cannot specify output name for list of files. This will only work with a single file (--in).");
    }
    if (!PVCalc && !GHCalc && normName == ""){
      _EXCEPTIONT("Need to either specify calculations for PV (--pv) or Z500 (--z500) or provide a name for the output variable (--normname).");
    }

   int nFiles,avgTime,nTime,nLat,nLon;
   if (const_thresh){
     if (PVCalc){
       anomVal = 1.2*std::pow(10,-6);
     }else if (GHCalc){
       anomVal = 170.;
     }else if (anomVal == 0.){
       _EXCEPTIONT("Anomaly threshold is currently set at 0. Specify a value using --threshval.");
     }
   }
   if (PVCalc){
     minThresh = 1.1*std::pow(10,-6);
     normName = "INT_ADVPV";
   }else if (GHCalc){
     minThresh = 100.;
     normName = "INT_ADZ";
   }

   if (const_thresh){
     outsuff = "norm_const.nc";
   }
   std::cout<<"Opening file list."<<std::endl;
   //Create list of input files
    std::vector<std::string> InputFiles;
  
    if (fileList != ""){
      GetInputFileList(fileList, InputFiles);
    }else{
      InputFiles.push_back(fileName);
    }
    nFiles = InputFiles.size();

    //Initialize a matrix. If constant, it will be filled with anom threshold values, else it will hold the 
    //threshold values

    std::cout<<"Opening reference file to load dimension variables."<<std::endl;
    int dim1= (endday-startday)+1;
    NcFile refFile(InputFiles[0].c_str());
    int dim2 = refFile.get_dim(latname.c_str())->size();
    int dim3 = refFile.get_dim(lonname.c_str())->size();
    refFile.close();

    DataArray3D <double> threshMat(dim1,dim2,dim3);

    if (!const_thresh){
      //Open threshold values file
      std::cout<<"Opening threshold file."<<std::endl;
      NcFile threshFile(threshName.c_str());
      if (!threshFile.is_valid()){
        _EXCEPTION1("Cannot open NetCDF file %s",threshName.c_str());
      }

      NcVar *threshVar = threshFile.get_var(threshVarName.c_str());
      if (threshVar == NULL){
        _EXCEPTION2("%s is missing variable %s",threshName.c_str(),threshVarName.c_str());
      }
      threshVar->set_cur(0,0,0);
      threshVar->get(&(threshMat[0][0][0]),dim1,dim2,dim3);
      threshFile.close();
    }else{
      for (int a=0; a<dim1; a++){
        for (int b=0; b<dim2; b++){
          for (int c=0; c<dim3; c++){
            threshMat[a][b][c] = anomVal;
          }
        }
      }
    }
    //Open var files

    for (int x=0; x<nFiles; x++){

      //Opening in read/write or read only depending on bool
      NcFile::FileMode readwrite;

      if (appendBool){
        readwrite=NcFile::Write;
        std::cout << "Appending variable to file "<< InputFiles[x].c_str() <<std::endl;
      }
      else{
        readwrite=NcFile::ReadOnly;
      }
      NcFile infile(InputFiles[x].c_str(),readwrite);
      if(!infile.is_valid()){
        _EXCEPTION1("Unable to open file \"%s\".",InputFiles[0].c_str());
      }
      std::cout<<"Opening file "<<InputFiles[x].c_str()<<std::endl;
      NcDim *tDim = infile.get_dim(tname.c_str());
      NcDim *latDim = infile.get_dim(latname.c_str());
      NcDim *lonDim = infile.get_dim(lonname.c_str());
      
      nTime = tDim->size();
      nLat = latDim->size();
      nLon = lonDim->size();
      NcVar *inTime = infile.get_var(tname.c_str());
      NcVar *inLat = infile.get_var(latname.c_str());
      NcVar *inLon = infile.get_var(lonname.c_str());
      NcVar *varData = infile.get_var(varName.c_str());

      NcAtt *attTime = inTime->get_att("units");
      if (attTime == NULL){
        _EXCEPTIONT("Time variable has no units attribute.");
      }

      std::string strTimeUnits = attTime->as_string(0);
      std::string strCalendar;
      NcAtt *attCal = inTime->get_att("calendar");
      if(attCal==NULL){
        strCalendar = "standard";
      }else{
        strCalendar = attCal->as_string(0);
      }


      DataArray1D<double> timeVals(nTime);
      inTime->set_cur((long) 0);
      inTime->get(&(timeVals[0]),nTime);

      //Create output file that corresponds to IPV data
  
      std::string strOutFile;
      if (outputName != ""){
        strOutFile = outputName;
      }else{
        pos = InputFiles[x].find(insuff);
        len = InputFiles[x].length();
        strOutFile = InputFiles[x].replace(pos,len,outsuff.c_str());
      }
    /*  std::string delim = ".";
      size_t pos,len;

      if (outputName != ""){
        strOutFile = outputName;
      }      
      else{
        pos = InputFiles[x].find(delim);
        len = InputFiles[x].length();
        strOutFile = InputFiles[x].replace(pos,len, "_norm.nc");
      }

      if (const_thresh){
        strOutFile = strOutFile.replace(strOutFile.end()-3,strOutFile.end(),"_const.nc");
      }*/
      if (appendBool){
        //Does this variable already exist within the file?
        NcVar *devIntOut = infile.get_var(normName.c_str());
        if (devIntOut == NULL){
          NcVar *devIntOut = infile.add_var(normName.c_str(),ncInt,tDim,latDim,lonDim);
        }
      //Create variables for Deviations
        if (PVCalc){
          calcNormalizedDevs(true,varData,devIntOut,inLat,inTime,strTimeUnits,strCalendar,threshMat,minThresh);

        }
        else{
          calcNormalizedDevs(false,varData,devIntOut,inLat,inTime,strTimeUnits,strCalendar,threshMat,minThresh);
        }
      }
      else{
        std::cout<<"Writing variables to file "<<strOutFile.c_str()<<std::endl;
        NcFile outfile(strOutFile.c_str(), NcFile::Replace, NULL,0,NcFile::Offset64Bits);
        int nOutTime = nTime;
        NcDim *tDimOut = outfile.add_dim(tname.c_str(),nOutTime);
        NcDim *latDimOut = outfile.add_dim(latname.c_str(),nLat);
        NcDim *lonDimOut = outfile.add_dim(lonname.c_str(),nLon);

        NcVar *tVarOut = outfile.add_var(tname.c_str(),ncDouble,tDimOut);
        CopyNcVarAttributes(inTime,tVarOut);      
        copy_dim_var(inTime,tVarOut);


        NcVar *latVarOut = outfile.add_var(latname.c_str(),ncDouble,latDimOut);
        copy_dim_var(inLat,latVarOut);
        NcVar *lonVarOut = outfile.add_var(lonname.c_str(),ncDouble,lonDimOut);
        copy_dim_var(inLon,lonVarOut);
        NcVar *devIntOut = outfile.add_var(normName.c_str(),ncInt,tDimOut,latDimOut,lonDimOut);
        //Create variables for Deviations

        if (PVCalc){
          calcNormalizedDevs(true,varData,devIntOut,inLat,tVarOut,strTimeUnits,strCalendar,threshMat,minThresh);
  
        }
        else{
          calcNormalizedDevs(false,varData,devIntOut,inLat,tVarOut,strTimeUnits,strCalendar,threshMat,minThresh);
        }
      }

    }
  }catch (Exception & e){
    std::cout<<e.ToString() <<std::endl;
  }

}
