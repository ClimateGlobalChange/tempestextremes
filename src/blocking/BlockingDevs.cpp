/////////////////////////////////////////////////////////
///
///           \file blockingDevs.cpp
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
#include "DataVector.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>


int main(int argc, char **argv){
//  NcError error(NcError::verbose_nonfatal);
  NcError error(NcError::silent_nonfatal);
  try{

    //list of input files for which to calculate deviations
    std::string fileList;
    //name of single input file
    std::string fileName;
    //name of single output file
    std::string outName;
    //file that holds averages
    std::string avgName;
    //Name of input variable
    std::string varName;
    //name of averaged input variable
    std::string avgVarName;
//    std::string threshName;
//    std::string threshVarName;
    std::string tname,latname,lonname;
    bool PVCalc;
    bool GHCalc;
//    bool const_thresh;


    BeginCommandLine()
      CommandLineString(fileName, "in","");
      CommandLineString(outName, "out","");
      CommandLineString(fileList, "inlist", "");
      CommandLineString(varName,"varname","");
      CommandLineString(avgName, "avg", "");
      CommandLineString(avgVarName, "avgname","");
//      CommandLineString(threshName,"thresh","");
//      CommandLineString(threshVarName,"threshname","");
      CommandLineBool(PVCalc,"pv");
      CommandLineBool(GHCalc,"gh");
//      CommandLineBool(const_thresh,"const");
      CommandLineString(tname,"tname","time");
      CommandLineString(latname,"latname","lat");
      CommandLineString(lonname,"lonname","lon");
      ParseCommandLine(argc, argv);
    EndCommandLine(argv);
    AnnounceBanner();

    if ((!PVCalc) && (!GHCalc)){
      _EXCEPTIONT("Need to specify either PV (--pv) or GH (--gh) calculations.");
    }
    if (fileName == "" && fileList == ""){
      _EXCEPTIONT("Need to specify either input file (--in) or file list (--inlist).");
    }
    if (fileName != "" && fileList != ""){
      _EXCEPTIONT("Cannot specify both file name (--in) and file list (--inlist).");
    }
    if (outName != "" && fileList != ""){
      _EXCEPTIONT("Currently cannot specify output name for list of files. This will only work with a single file (--in).");
    }

   int nFiles,avgTime,nTime,nLat,nLon;
   double anomVal;
   if (PVCalc){
     anomVal = 1.3*std::pow(10,-6);
   }else if (GHCalc){
     anomVal = 170.;
   }
   //Create list of input files
    std::vector<std::string> InputFiles;

    if (fileList != ""){
      GetInputFileList(fileList, InputFiles);
    }
    else{
      InputFiles.push_back(fileName);
    }
    nFiles = InputFiles.size();

    //Open averages file
    NcFile avgFile(avgName.c_str());
	if (!avgFile.is_valid()) {
		_EXCEPTION1("Cannot open NetCDF file \"%s\"", avgName.c_str());
	}

    //time data (average)
	NcDim *dimTime = avgFile.get_dim(tname.c_str());
	if (dimTime == NULL) {
		_EXCEPTION1("\"%s\" is missing dimension \"time\"", avgName.c_str());
	}
    avgTime = dimTime->size();
    NcVar *avgTimeVals = avgFile.get_var(tname.c_str());
   	if (avgTimeVals == NULL) {
		_EXCEPTION1("\"%s\" is missing variable \"time\"", avgName.c_str());
	}
    
    //averaged var
    NcVar *AvarData = avgFile.get_var(avgVarName.c_str());
   	if (AvarData == NULL) {
		_EXCEPTION2("\"%s\" is missing variable \"%s\"", avgName.c_str(), avgVarName.c_str());
	}
    int dim1 = AvarData->get_dim(0)->size();
    int dim2 = AvarData->get_dim(1)->size();
    int dim3 = AvarData->get_dim(2)->size();    
/*    //Initialize a matrix. If constant, it will be filled with anom values, else it will hold the 
    //threshold values

    DataMatrix3D <double> threshMat(dim1,dim2,dim3);

    if (!const_thresh){
      //Open threshold values file

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
*/
    //Open var files

    for (int x=0; x<nFiles; x++){

      NcFile infile(InputFiles[x].c_str());
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
   //   std::cout<<"Time units: "<< strTimeUnits<<" Calendar: "<<strCalendar<<std::endl;



      //Get the day index of the first time step
      DataVector<double> timeVals(nTime);
      inTime->set_cur((long) 0);
      inTime->get(&(timeVals[0]),nTime);

      int dateYear; 
      int dateMonth;
      int dateDay;
      int dateHour;
      //std::cout<<"time units and calendar: "<<strTimeUnits<<","<<strCalendar<<std::endl;
      //std::cout<<"first value of time variable:"<<timeVals[0]<<std::endl;
      ParseTimeDouble(strTimeUnits, strCalendar, timeVals[0], dateYear,\
        dateMonth, dateDay, dateHour);
      //std::cout<<"D/M/Y:"<<dateDay<<"/"<<dateMonth<<"/"<<dateYear<<std::endl;
      int day = DayInYear(dateMonth,dateDay);
      //std::cout<<"For month "<<dateMonth<<" and day "<<dateDay<<" day is "<<day<<std::endl;
      int startIndex = day-1;

    //  int nSteps = int(1.0/(timeVals[1]-timeVals[0]));
      double tRes;
      if ((strTimeUnits.length() >= 11) && \
        (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)){
        tRes = timeVals[1]-timeVals[0];
      }
      else if((strTimeUnits.length() >= 12) && \
      (strncmp(strTimeUnits.c_str(), "hours since ",12)==0)) {
        tRes = (timeVals[1]-timeVals[0])/24.0;
      }else if (strTimeUnits.length() >= 14 && \
        (strncmp(strTimeUnits.c_str(),"minutes since ",14)== 0) ){
        tRes = (timeVals[1]-timeVals[0])/(24. * 60.);
      }else{
       _EXCEPTIONT("Cannot determine time resolution (unknown time units).");
      }
      int nSteps = 1/tRes;

    //check if the file contains leap days 
      bool leap = false;
      int leapYear = 0;
      int leapMonth = 0;
      int leapDay = 0;
      int leapHour = 0;

      int nLeapSteps = 0;
      if (strCalendar!="noleap"){
        for (int t=0; t<nTime; t++){
          ParseTimeDouble(strTimeUnits, strCalendar, timeVals[t], leapYear,\
            leapMonth, leapDay, leapHour);
          if ((leapMonth==2 && leapDay == 29)){
            leap = true;
            nLeapSteps +=1;
          }     
        }
      }
      

      //Create output file that corresponds to IPV data
      std::string strOutFile;
      if (outName != ""){
        strOutFile = outName;
      }else{
        strOutFile = InputFiles[x].replace(InputFiles[x].end()-3,\
          InputFiles[x].end(), "_devs.nc");
      }
      std::cout<<"Writing variables to file "<<strOutFile.c_str()<<std::endl;
      NcFile outfile(strOutFile.c_str(), NcFile::Replace, NULL,0,NcFile::Offset64Bits);
      int nOutTime;
      nOutTime = nTime-nLeapSteps;
      NcDim *tDimOut = outfile.add_dim(tname.c_str(),nOutTime);
      NcDim *latDimOut = outfile.add_dim(latname.c_str(),nLat);
      NcDim *lonDimOut = outfile.add_dim(lonname.c_str(),nLon);

      NcVar *tVarOut = outfile.add_var(tname.c_str(),ncDouble,tDimOut);
      CopyNcVarAttributes(inTime,tVarOut);      

      NcVar *latVarOut = outfile.add_var(latname.c_str(),ncDouble,latDimOut);
      copy_dim_var(inLat,latVarOut);
      NcVar *lonVarOut = outfile.add_var(lonname.c_str(),ncDouble,lonDimOut);
      copy_dim_var(inLon,lonVarOut);

      //Create variables for Deviations

      if (PVCalc){
        NcVar *devOut = outfile.add_var("DIPV",ncDouble,tDimOut,latDimOut,lonDimOut);
        NcVar *aDevOut = outfile.add_var("ADIPV",ncDouble,tDimOut,latDimOut,lonDimOut);
  //      NcVar *devIntOut = outfile.add_var("INT_ADIPV",ncInt,tDimOut,latDimOut,lonDimOut);

        calcDevs(leap,true, startIndex, tRes, strTimeUnits, strCalendar, varData, devOut, aDevOut, AvarData, inTime,\
        avgTimeVals, inLat, tVarOut);
      }
      else if (GHCalc){
        NcVar *devOut = outfile.add_var("DGH",ncDouble,tDimOut,latDimOut,lonDimOut);
        NcVar *aDevOut = outfile.add_var("ADGH",ncDouble,tDimOut,latDimOut,lonDimOut);
//        NcVar *devIntOut = outfile.add_var("INT_ADGH",ncInt,tDimOut,latDimOut,lonDimOut);
//        NcVar *stdDevOut = outfile.add_var("STD_DEV",ncDouble,latDimOut,lonDimOut);
        calcDevs(leap,false, startIndex, tRes, strTimeUnits, strCalendar, varData, devOut,aDevOut,AvarData,inTime,\
          avgTimeVals,inLat,tVarOut);
        std::cout<<"Finished writing to file "<<strOutFile.c_str()<<std::endl;
      }
      else{
        _EXCEPTIONT("Invalid variable specified!");
      }
    }
  }catch (Exception & e){
    std::cout<<e.ToString() <<std::endl;
  }

}
