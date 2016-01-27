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

#include "blockingUtilities.h"
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

void calcDevs(bool leap,
              int startAvgIndex,
              NcVar *inIPV, 
              NcVar *outDev, 
              NcVar *outADev,
              NcVar *outPosIntDev,
              NcVar *avgIPV, 
              NcVar *inTime, 
              NcVar *avgTime,
              NcVar *lat,
              NcVar *outTime,
              double PVAnom){

  int nTime,nLat,nLon,nSteps,avgDay,nOutTime;
  double tRes;

  nTime = inIPV->get_dim(0)->size();
  nLat = inIPV->get_dim(1)->size();
  nLon = inIPV->get_dim(2)->size();

//keep sign consistent!
  if (PVAnom<0){
    PVAnom = -PVAnom;
  }

//input IPV
  DataMatrix3D<double> IPVMat(nTime,nLat,nLon);
  inIPV->set_cur(0,0,0);
  inIPV->get(&(IPVMat[0][0][0]),nTime,nLat,nLon);

  DataVector<double> timeVec(nTime);
  inTime->set_cur((long) 0);
  inTime->get(&(timeVec[0]),nTime);

  std::string strTimeUnits = inTime->get_att("units")->as_string(0);
  std::string strCalendar = inTime->get_att("calendar")->as_string(0);

  if ((strTimeUnits.length() >= 11) && \
    (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)){
    tRes = timeVec[1]-timeVec[0];
  }
  else{
    tRes = (timeVec[1]-timeVec[0])/24.0;
  }
  nSteps = 1/tRes;
 
//avg IPV
  avgDay = avgIPV->get_dim(0)->size();
  DataMatrix3D<double> avgMat(avgDay,nLat,nLon);
  avgIPV->set_cur(0,0,0);
  avgIPV->get(&(avgMat[0][0][0]),avgDay,nLat,nLon);

  DataVector<int> avgDayVec(avgDay);
  avgTime->set_cur((long) 0);
  avgTime->get(&(avgDayVec[0]),avgDay);


//Matrix for output data
//Eliminate one day if contains Feb 29
  if (leap){
    nOutTime = nTime-nSteps;
    std::cout<<"# steps previously "<<nTime<<" but now "<<nOutTime<<std::endl;
  }
  else{
    nOutTime = nTime;
  }

  DataMatrix3D<double> devMat(nOutTime,nLat,nLon);
  DataMatrix3D<double> aDevMat(nOutTime,nLat,nLon);

//Number of days in IPV
  int nDays = nTime*tRes;
  std::cout<<"There are "<<nDays<<" days in file."<<std::endl;


//Deal with skipped days          
  int d=0;
  DataVector<double> newTime(nOutTime);

  int leapYear=0;
  int leapMonth=0;
  int leapDay=0;
  int leapHour=0;

  for (int t=0; t<nTime; t++){
    if (leap){
      ParseTimeDouble(strTimeUnits, strCalendar, timeVec[t], leapYear,\
        leapMonth, leapDay, leapHour);
      if (leapMonth==2 && leapDay == 29){
       /* while (leapMonth ==2 && leapDay == 29){
          std::cout<<"Leap day! Skipping this step."<<std::endl;
          t++;
          ParseTimeDouble(strTimeUnits,strCalendar, timeVec[t],leapYear,\
            leapMonth,leapDay,leapHour);
        }*/
        std::cout<<"Leap day! Skipping day."<<std::endl;
        t+=nSteps;
      }      
    }
    if (t>=nTime){
      std::cout<<"t equals or exceeds nTime. Leaving loop."<<std::endl;
      break;
    }
    int nDayIncrease = d/nSteps;
    std::cout<<"t is "<<t<<" and current number of days past start is "<<nDayIncrease<<std::endl;
    int currAvgIndex = startAvgIndex + nDayIncrease;
    if (currAvgIndex>364){
      currAvgIndex-=365;
    }
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        devMat[t][a][b] = IPVMat[t][a][b]-avgMat[currAvgIndex][a][b];
      }
    }
    newTime[d] = timeVec[t];
    std::cout<<"Day is "<<currAvgIndex<<", d is "<<d<<" and t is "<<t<<std::endl;
    d++;
  }
  outTime->set_cur((long) 0);
  outTime->put(&(newTime[0]),nOutTime);

  std::cout<<"About to implement smoothing."<<std::endl;
  double div = (double) 2*nSteps;
  double invDiv = 1.0/div;

  //implement 2-day smoothing
  for (int t=0; t<nOutTime; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        if (t<2*nSteps){
          aDevMat[t][a][b] = devMat[t][a][b];
        }
        else{
          for (int n=0; n<2*nSteps; n++){
            aDevMat[t][a][b]+=devMat[t-n][a][b];
          }
          aDevMat[t][a][b] = aDevMat[t][a][b]*invDiv; 
        }
      }
    }
  }

  std::cout<<"Finished smoothing."<<std::endl;
  outDev->set_cur(0,0,0);
  outDev->put(&(devMat[0][0][0]),nOutTime,nLat,nLon);
  std::cout<<"Wrote devs to file."<<std::endl;

  outADev->set_cur(0,0,0);
  outADev->put(&(aDevMat[0][0][0]),nOutTime,nLat,nLon);
  std::cout<<"Wrote smoothed devs to file."<<std::endl;
//Divide matrix by PV anomaly value 
//We are looking for negative anomalies in NH and positive anomalies in SH
  DataVector<double> latVec(nLat);
  lat->set_cur((long) 0);
  lat->get(&(latVec[0]),nLat);

  DataMatrix3D<int> posIntDevs(nOutTime,nLat,nLon);
  double invAnom = 1.0/PVAnom;
  double divDev,pos,neg;
  for (int t=0; t<nOutTime; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        divDev = aDevMat[t][a][b]*invAnom;
        //SH: positive anomalies
        if (latVec[a]<0){
          pos = (divDev+std::fabs(divDev))*0.5;
          posIntDevs[t][a][b] = int(pos);
        }
        //NH: negative anomalies
        else if (latVec[a]>=0){
          neg = (divDev-std::fabs(divDev))*0.5;
          posIntDevs[t][a][b] = -int(neg);
        }
      }
    }
  }
  outPosIntDev->set_cur(0,0,0);
  outPosIntDev->put(&(posIntDevs[0][0][0]),nOutTime,nLat,nLon);
  std::cout<<"Wrote integer values to file."<<std::endl;
}





int main(int argc, char **argv){
  NcError error(NcError::verbose_nonfatal);

  try{

    //list of input files for which to calculate deviations
    std::string fileList;
    //file that holds averages
    std::string avgName;

    BeginCommandLine()
      CommandLineString(fileList, "inlist", "");
      CommandLineString(avgName, "avg", "");
      ParseCommandLine(argc, argv);
    EndCommandLine(argv);
    AnnounceBanner();

   int nFiles,avgTime,nTime,nLat,nLon;
   double anomVal = 1.3*std::pow(10,-6);

   //Create list of input files
    std::vector<std::string> InputFiles;
    GetInputFileList(fileList, InputFiles);
    nFiles = InputFiles.size();

    //Open averages file
    NcFile avgFile(avgName.c_str());

    //time data (average)
    avgTime = avgFile.get_dim("time")->size();
    NcVar *avgTimeVals = avgFile.get_var("time");
    
    //averaged IPV
    NcVar *AIPVData = avgFile.get_var("AIPV");

    //Open IPV files

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
      NcVar *IPVData = infile.get_var("IPV");

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
      NcVar *devOut = outfile.add_var("DIPV",ncDouble,tDimOut,latDimOut,lonDimOut);
      NcVar *aDevOut = outfile.add_var("ADIPV",ncDouble,tDimOut,latDimOut,lonDimOut);
      NcVar *devIntOut = outfile.add_var("INT_ADIPV",ncInt,tDimOut,latDimOut,lonDimOut);

      calcDevs(leap, startIndex, IPVData, devOut, aDevOut, devIntOut, AIPVData, inTime,\
        avgTimeVals, inLat, tVarOut, anomVal);
 
      infile.close();
      outfile.close();
    }

  }catch (Exception & e){
    std::cout<<e.ToString() <<std::endl;
  }

}