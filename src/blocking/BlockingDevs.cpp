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
  NcError error(NcError::verbose_nonfatal);

  try{

    //list of input files for which to calculate deviations
    std::string fileList;
    //file that holds averages
    std::string avgName;
    std::string varName;
    std::string avgVarName;

    BeginCommandLine()
      CommandLineString(fileList, "inlist", "");
      CommandLineString(varName,"varname","");
      CommandLineString(avgName, "avg", "");
      CommandLineString(avgVarName, "avgname","");
      ParseCommandLine(argc, argv);
    EndCommandLine(argv);
    AnnounceBanner();

   int nFiles,avgTime,nTime,nLat,nLon;
   double anomVal = 1.3*std::pow(10,-6);
   double GHVal = 200;
   //Create list of input files
    std::vector<std::string> InputFiles;
    GetInputFileList(fileList, InputFiles);
    nFiles = InputFiles.size();

    //Open averages file
    NcFile avgFile(avgName.c_str());

    //time data (average)
    avgTime = avgFile.get_dim("time")->size();
    NcVar *avgTimeVals = avgFile.get_var("time");
    
    //averaged var
    NcVar *AvarData = avgFile.get_var(avgVarName.c_str());

    //Open var files

    for (int x=0; x<nFiles; x++){

      NcFile infile(InputFiles[x].c_str());
      if(!infile.is_valid()){
        _EXCEPTION1("Unable to open file \"%s\".",InputFiles[0].c_str());
      }
      NcDim *tDim = infile.get_dim("time");
      NcDim *latDim = infile.get_dim("lat");
      NcDim *lonDim = infile.get_dim("lon");
      
      nTime = tDim->size();
      nLat = latDim->size();
      nLon = lonDim->size();

      NcVar *inTime = infile.get_var("time");
      NcVar *inLat = infile.get_var("lat");
      NcVar *inLon = infile.get_var("lon");
      NcVar *varData = infile.get_var(varName.c_str());

      NcAtt *attTime = inTime->get_att("units");
      if (attTime == NULL){
        _EXCEPTIONT("Time variable has no units attribute.");
      }

      std::string strTimeUnits = attTime->as_string(0);

      NcAtt *attCal = inTime->get_att("calendar");
      if(attCal==NULL){
        _EXCEPTIONT("Time variable has no calendar attribute.");
      }
      std::string strCalendar = attCal->as_string(0);

      std::cout<<"Time units: "<< strTimeUnits<<" Calendar: "<<strCalendar<<std::endl;



      //check if file is leap year file
      bool leap = false;
      DataVector<double> timeVals(nTime);
      inTime->set_cur((long) 0);
      inTime->get(&(timeVals[0]),nTime);

      int dateYear; 
      int dateMonth;
      int dateDay;
      int dateHour;

      ParseTimeDouble(strTimeUnits, strCalendar, timeVals[0], dateYear,\
        dateMonth, dateDay, dateHour);

      int day = DayInYear(dateMonth,dateDay);
      int startIndex = day-1;

    //  int nSteps = int(1.0/(timeVals[1]-timeVals[0]));
      double tRes;
      if ((strTimeUnits.length() >= 11) && \
        (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)){
        tRes = timeVals[1]-timeVals[0];
      }
      else{
        tRes = (timeVals[1]-timeVals[0])/24.0;
      }
      int nSteps = 1/tRes;
     
      int leapYear = 0;
      int leapMonth = 0;
      int leapDay = 0;
      int leapHour = 0;

      if (strCalendar!="noleap" && dateMonth <=2){
        ParseTimeDouble(strTimeUnits, strCalendar, timeVals[nTime-1], leapYear,\
          leapMonth, leapDay, leapHour);
        if ((leapMonth==2 && leapDay == 29) || (dateMonth ==2 && leapMonth == 3)){
          leap = true;
        }     

      }
      

      //Create output file that corresponds to IPV data
      std::string strOutFile = InputFiles[x].replace(InputFiles[x].end()-3,\
        InputFiles[x].end(), "_devs.nc");

      NcFile outfile(strOutFile.c_str(), NcFile::Replace, NULL,0,NcFile::Offset64Bits);
      int nOutTime;
      if (leap){
        nOutTime = nTime-nSteps;
      }
      else{
        nOutTime = nTime;
      }
      NcDim *tDimOut = outfile.add_dim("time",nOutTime);
      NcDim *latDimOut = outfile.add_dim("lat",nLat);
      NcDim *lonDimOut = outfile.add_dim("lon",nLon);

      NcVar *tVarOut = outfile.add_var("time",ncDouble,tDimOut);
      CopyNcVarAttributes(inTime,tVarOut);      

      NcVar *latVarOut = outfile.add_var("lat",ncDouble,latDimOut);
      copy_dim_var(inLat,latVarOut);
      NcVar *lonVarOut = outfile.add_var("lon",ncDouble,lonDimOut);
      copy_dim_var(inLon,lonVarOut);

      //Create variables for Deviations

      if (varName == "IPV"){
        NcVar *devOut = outfile.add_var("DIPV",ncDouble,tDimOut,latDimOut,lonDimOut);
        NcVar *aDevOut = outfile.add_var("ADIPV",ncDouble,tDimOut,latDimOut,lonDimOut);
        NcVar *devIntOut = outfile.add_var("INT_ADIPV",ncInt,tDimOut,latDimOut,lonDimOut);

        calcDevsPV(leap, startIndex, varData, devOut, aDevOut, devIntOut, AvarData, inTime,\
        avgTimeVals, inLat, tVarOut, anomVal);
      }
      else if (varName == "GH"){
        NcVar *devOut = outfile.add_var("DGH",ncDouble,tDimOut,latDimOut,lonDimOut);
        NcVar *aDevOut = outfile.add_var("ADGH",ncDouble,tDimOut,latDimOut,lonDimOut);
        NcVar *devIntOut = outfile.add_var("INT_ADGH",ncInt,tDimOut,latDimOut,lonDimOut);
        NcVar *stdDevOut = outfile.add_var("STD_DEV",ncDouble,latDimOut,lonDimOut);
        calcDevsGH(leap, startIndex, varData, devOut,aDevOut,devIntOut,AvarData,inTime,avgTimeVals,inLat,tVarOut,stdDevOut);
      }
      else{
        _EXCEPTIONT("Invalid variable specified!");
      }
    }
  }catch (Exception & e){
    std::cout<<e.ToString() <<std::endl;
  }

}
