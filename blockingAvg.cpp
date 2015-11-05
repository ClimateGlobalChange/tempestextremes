/////////////////////////////////////
///     \file blockingAvg.cpp
///     \author Marielle Pinheiro
///     \version March 26, 2015


/*This code fills an array with 31 days (31*nsteps) of PV
data 

*/
#include "blockingUtilities.h"
#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "TimeObj.h"

#include "DataVector.h"
#include "DataMatrix.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <sstream>
int main(int argc, char **argv){
  NcError error(NcError::verbose_nonfatal);

  try{
  std::string fileList;
  std::string strfile_out;

  BeginCommandLine()
    CommandLineString(fileList, "inlist", "");
    CommandLineString(strfile_out, "out", "");
    ParseCommandLine(argc, argv);

  EndCommandLine(argv)
  AnnounceBanner();

  //Create list of input files
  std::vector<std::string> InputFiles;
  GetInputFileList(fileList, InputFiles);
  int nFiles = InputFiles.size();

  //Open first file 
  NcFile infile(InputFiles[0].c_str());
  int nTime = infile.get_dim("time")->size();
  int nLat = infile.get_dim("lat")->size();
  int nLon = infile.get_dim("lon")->size();

  //IPV variable and data matrix
  NcVar *inPV = infile.get_var("IPV");

  DataMatrix3D<double> IPVData(nTime, nLat, nLon);
  inPV->set_cur(0,0,0);
  inPV->get(&(IPVData[0][0][0]),nTime,nLat,nLon);

  //time variable
  NcVar *timeVal = infile.get_var("time");

  DataVector<double> timeVec(nTime);
  timeVal->set_cur((long) 0);
  timeVal->get(&(timeVec[0]),nTime);

  NcAtt *attTime = timeVal->get_att("units");
  if (attTime == NULL){
    _EXCEPTIONT("Time variable has no units attribute.");
  }

  std::string strTimeUnits = attTime->as_string(0);

  NcAtt *attCal = timeVal->get_att("calendar");
  if(attCal==NULL){
    _EXCEPTIONT("Time variable has no calendar attribute.");
  }
  std::string strCalendar = attCal->as_string(0);
  if (strncmp(strCalendar.c_str(), "gregorian",9)==0){
    std::cout<< "Changing calendar type to standard."<<std::endl;
    strCalendar = "standard";
  }

  double tRes;
  if ((strTimeUnits.length() >= 11) && \
    (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)){
    tRes = timeVec[1]-timeVec[0];
  }
  else {
    tRes = (timeVec[1]-timeVec[0])/24.0;
  }

//Parse time units of first file to determine start date 
  int
 yearLen = 365;
  int dateYear=0;
  int dateMonth=0;
  int dateDay=0;
  int dateHour=0;

  //Length of 31 days axis
  int arrLen = int(31.0/tRes);

  //3D matrix to store averaged values
  DataMatrix3D<double> avgStoreVals(yearLen,nLat,nLon);
  DataMatrix3D<double> avgCounts(yearLen,nLat,nLon);

  //Start date of first file
  ParseTimeDouble(strTimeUnits, strCalendar, timeVec[0], dateYear,\
    dateMonth, dateDay, dateHour);

  int day = DayInYear(dateMonth,dateDay);
  int dateIndex = day + 15;

  int leapYear=0;
  int leapMonth=0;
  int leapDay=0;
  int leapHour=0;

  bool leap = false;
  if (strCalendar!="noleap" && dateMonth<=2){
    //Check whether file contains a Feb 29

    ParseTimeDouble(strTimeUnits, strCalendar, timeVec[nTime-1], leapYear,\
      leapMonth, leapDay, leapHour);

    if ((leapMonth==2 && leapDay==29) || (dateMonth==2&&leapMonth==3)){
      //Check when parsing the indices
      std::cout<<"Might contain leap day. Will check."<<std::endl;
      leap = true;
    }
  }

  //Number of time steps per day
  int nSteps = 1/tRes;
  std::cout<<"tRes is "<<tRes<<" and nSteps is "<<nSteps<<std::endl;
  double endTime = timeVec[nTime-1];
  int currArrIndex = 0;

  //3D 31-day array
  DataMatrix3D<double> currFillData(arrLen,nLat,nLon);

  int tEnd;
  //File counter
  int x=0;
  //Start and end for time dimension
  int tStart = 0;
  if (nTime>arrLen){
    tEnd = arrLen;
  }
  else{
    tEnd = nTime;
  }

  //First while loop: open files and fill until 31 days array full
  while (currArrIndex<arrLen){
    for (int t=tStart; t<tEnd; t++){
      if (leap==true){
      //Check if time is a leap year date
        ParseTimeDouble(strTimeUnits, strCalendar, timeVec[t], leapYear,\
          leapMonth,leapDay,leapHour);
        if (leapMonth==2 && leapDay == 29){
          while (leapMonth ==2 && leapDay ==29){
            std::cout<<"Leap day! Skipping this step."<<std::endl;
            t++;
            ParseTimeDouble(strTimeUnits, strCalendar, timeVec[t], leapYear,\
              leapMonth, leapDay, leapHour);
          }
        }
      }
      //Fill the 31 day array with PV data
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          currFillData[currArrIndex][a][b] = IPVData[t][a][b];
        }
      }
      currArrIndex++;
    }
    //If the 31 days aren't yet filled, keep going
    if (currArrIndex<arrLen){
      infile.close();
      //Open new file and increment file counter
      x+=1;
      NcFile infile(InputFiles[x].c_str());
      nTime = infile.get_dim("time")->size();
      nLat = infile.get_dim("lat")->size();
      nLon = infile.get_dim("lon")->size();

      //IPV variable
      NcVar *inPV = infile.get_var("IPV");
      IPVData.Initialize(nTime,nLat,nLon);

      inPV->set_cur(0,0,0);
      inPV->get(&(IPVData[0][0][0]),nTime,nLat,nLon);

      //Time variable
      NcVar *timeVal = infile.get_var("time");
      timeVec.Initialize(nTime);
      timeVal->set_cur((long) 0);
      timeVal->get(&(timeVec[0]),nTime);
      
      //check to make sure that there isn't a file missing in the list!
      //otherwise will mess up averages
      double contCheck;
      
      if ((strTimeUnits.length() >= 11) && \
        (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)){
        contCheck = std::fabs(timeVec[0]-endTime);
      }
      else{
        contCheck = std::fabs(timeVec[0]-endTime)/24.0;
      }

      if (contCheck>tRes){
        _EXCEPTIONT("New file is not continuous with previous file."); 
      }

      //reset leap year check
      leap = false;

      //check for leap year days in this file
      ParseTimeDouble(strTimeUnits, strCalendar, timeVec[0], dateYear,\
        dateMonth, dateDay, dateHour);

      if (strCalendar!="noleap" && dateMonth<=2){
      //Check whether file contains a Feb 29

        ParseTimeDouble(strTimeUnits, strCalendar, timeVec[nTime-1], leapYear,\
          leapMonth, leapDay, leapHour);

        if ((leapMonth==2 && leapDay==29) || (dateMonth==2 && leapMonth==3)){
        //Check when parsing the indices
          std::cout<<"May contain leap day. Will check."<<std::endl;
          leap = true;
        }
      }
      //reset ending time of current file for next continuity check
      endTime = timeVec[nTime-1];

      //check that tEnd doesn't exceed 31 days
     
      int tCheck = currArrIndex + nTime;

      if (tCheck > arrLen){
        tEnd = arrLen-currArrIndex;
      }
      else{
        tEnd = nTime;
      }
      std::cout<<"tEnd is "<<tEnd<<std::endl;
    }
    
  }
  //Fill yearly array with sum of 31 days
  for (int t=0; t<arrLen; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        avgStoreVals[dateIndex][a][b] += currFillData[t][a][b];
      }
    }
  }
 //increase count by 1
  for (int a=0; a<nLat; a++){
    for (int b=0; b<nLon; b++){
      avgCounts[dateIndex][a][b] +=1.0;
    }
  }
  
  dateIndex+=1;
  currArrIndex = 0;

//Check if new file needs to be opened before entering while loop

  if(tEnd>=nTime){
    std::cout<<"tEnd is "<<tEnd<<". Reached end of previous file."<<std::endl;
    infile.close();
    std::cout<<"Closed "<<InputFiles[x]<<std::endl;
    x+=1;
    NcFile infile(InputFiles[x].c_str());
    nTime = infile.get_dim("time")->size();
    nLat = infile.get_dim("lat")->size();
    nLon = infile.get_dim("lon")->size();

    //IPV variable
    NcVar *inPV = infile.get_var("IPV");
    IPVData.Initialize(nTime,nLat,nLon);

    inPV->set_cur(0,0,0);
    inPV->get(&(IPVData[0][0][0]),nTime,nLat,nLon);

    //Time variable
    NcVar *timeVal = infile.get_var("time");
    timeVec.Initialize(nTime);

    timeVal->set_cur((long) 0);
    timeVal->get(&(timeVec[0]),nTime);

    leap = false;
 
    ParseTimeDouble(strTimeUnits, strCalendar, timeVec[0], dateYear,\
      dateMonth, dateDay, dateHour);

    if (strCalendar!="noleap" && dateMonth<=2){
      //Check whether file contains a Feb 29
      std::cout<<"Checking leap year status."<<std::endl;
      ParseTimeDouble(strTimeUnits, strCalendar, timeVec[nTime-1], leapYear,\
        leapMonth, leapDay, leapHour);
      if ((leapMonth==2 && leapDay==29) || (dateMonth==2 && leapMonth==3)){
      //Check when parsing the indices
        std::cout<<"May contain leap day. Will check."<<std::endl;
        leap = true;
      }
    }
    tStart = 0;
    tEnd = tStart + nSteps;
  }

  else if (tEnd<nTime){
    std::cout<<"Still have data left on current file. Will continue on."<<std::endl;
    tStart = tEnd;
    tEnd = tStart + nSteps;
  }

  //entering while loop.
  std::cout<<"Before entering: start, end is "<<tStart<<", "<<tEnd<<std::endl;
  bool newFile = false;

  while (x<nFiles){
    std::cout<<"7: Inside while loop: file currently "<<InputFiles[x]<<" and nTime is "\
      <<nTime<<"; tStart/end is "<<tStart<<" and "<<tEnd<<std::endl;
     //1: check if current day is a leap day
    if (leap==true){
      std::cout<<"Checking date for leap day."<<std::endl;
      ParseTimeDouble(strTimeUnits, strCalendar, timeVec[tStart], leapYear,\
        leapMonth, leapDay, leapHour);
      if (leapMonth==2 && leapDay ==29){
        tStart = tEnd;
        tEnd = tStart + nSteps;
        std::cout<<"This day is a leap day. Resetting tStart/tEnd to "\
          <<tStart<<" and "<<tEnd<<std::endl;
      } 
      if (tEnd>nTime){
        newFile = true;
        std::cout<<"For leap day file, reached EOF."<<std::endl;
      }
    }
    //2: if new file needs to be opened, open it
    if(newFile==true){
      std::cout<<"Reached end of file."<<std::endl;
      infile.close();
      std::cout<<"8: Closed "<<InputFiles[x]<<std::endl;
      x+=1;
      if (x<(nFiles-1)){

        NcFile infile(InputFiles[x].c_str());
        newFile = false;
        std::cout<<"x is currently "<<x<<", opening file "<<InputFiles[x]<<std::endl;
        nTime = infile.get_dim("time")->size();
        nLat = infile.get_dim("lat")->size();
        nLon = infile.get_dim("lon")->size();

      //IPV variable
        NcVar *inPV = infile.get_var("IPV");
        IPVData.Initialize(nTime, nLat, nLon);

        inPV->set_cur(0,0,0);
        inPV->get(&(IPVData[0][0][0]),nTime,nLat,nLon);

      //Time variable
        NcVar *timeVal = infile.get_var("time");
        timeVec.Initialize(nTime);
        timeVal->set_cur((long) 0);
        timeVal->get(&(timeVec[0]),nTime);

        leap = false;

        ParseTimeDouble(strTimeUnits, strCalendar, timeVec[0], dateYear,\
          dateMonth, dateDay, dateHour);

        if (strCalendar!="noleap" && dateMonth<=2){
        //Check whether file contains a Feb 29
          std::cout<<"Checking leap year status."<<std::endl;
          ParseTimeDouble(strTimeUnits, strCalendar, timeVec[nTime-1], leapYear,\
          leapMonth, leapDay, leapHour);

          if ((leapMonth==2 && leapDay==29) || (dateMonth==2 && leapMonth==3)){
        //Check when parsing the indices
            std::cout<<"May contain leap day. Will check."<<std::endl;
            leap = true;
          }
        }
        tStart = 0;
        tEnd = tStart + nSteps;
      }
    }    
    //3: fill array for next value
    if (newFile==false){
  //    std::cout<<"DEBUG for fill: tStart, tEnd are "<<tStart<<", "<<tEnd\
        <<"and nTime is "<<nTime<<". Now parsing dates for those times."<<std::endl;


      for (int t=tStart; t<tEnd; t++){
        for (int a=0; a<nLat; a++){
          for (int b=0; b<nLon; b++){
            currFillData[currArrIndex][a][b] = IPVData[t][a][b];
          }
        }
        currArrIndex+=1;
      }
      //check for periodic boundary condition for 31 day array
      if (currArrIndex>=arrLen){
        std::cout<<"currArrIndex is "<<currArrIndex<<" and periodic boundary: 31 day length met or exceeded."<<std::endl;
        currArrIndex-=arrLen;
        std::cout<<"currArrIndex is now "<<currArrIndex<<std::endl;
      }
   //Fill date with sum of new array values
      for (int t=0; t<arrLen; t++){
        for (int a=0; a<nLat; a++){
          for (int b=0; b<nLon; b++){
            avgStoreVals[dateIndex][a][b] += currFillData[t][a][b];
          }
        }
      }
      //increase count
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          avgCounts[dateIndex][a][b] += 1.0;
          if (a==50 && b==50){
            std::cout<<" date index is "<<dateIndex<<" and count is "<<avgCounts[dateIndex][a][b]<<std::endl;
          }
        }
      }
      dateIndex+=1;
      //periodic boundary condition for year
      if(dateIndex>=yearLen){
        dateIndex-=yearLen;
        std::cout<<"Periodic boundary: date index is "<<dateIndex<<" and end of year array reached. dateIndex is now "\
          <<dateIndex<<std::endl;
      }

      //
      if (tEnd<nTime){
        std::cout<<"start/end/ntime:"<<tStart<<"/"<<tEnd<<"/"<<nTime<<". Still have data left on current file. Will continue on."<<std::endl;
        tStart = tEnd;
        tEnd = tStart + nSteps;
      //  std::cout<<"start/end now "<<tStart <<" and "<< tEnd <<std::endl;
        //Check for leap year
        if (leap==true){
          std::cout<<"Checking date for leap day."<<std::endl;
          ParseTimeDouble(strTimeUnits, strCalendar, timeVec[tStart], leapYear,\
            leapMonth, leapDay, leapHour);
          if (leapMonth==2 && leapDay ==29){
            tStart = tEnd;
            tEnd = tStart + nSteps;
            std::cout<<"This day is a leap day. Resetting tStart/tEnd to "\
              <<tStart<<" and "<<tEnd<<std::endl;
          }
          //re-check tEnd
          if (tEnd>=nTime){
            std::cout<<"After re-check of tEnd (leap day), reached EOF"<<std::endl;
            newFile = true;
          } 
        }
      }
      else{
        newFile = true;
        std::cout<<"Else: Reached EOF. Will open new file in next iteration."<<std::endl;
      }
    }
  }
  
  //average all values
  for (int t=0; t<yearLen; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        if (a==50 && b==50){
          std::cout<<"DEBUG: for t="<<t<<", value is "<<avgStoreVals[t][a][b]\
            <<", count is "<<avgCounts[t][a][b]<<" and dividing by "<<avgCounts[t][a][b]*31.0*nSteps<<std::endl;
        }
        avgStoreVals[t][a][b] = avgStoreVals[t][a][b]/(avgCounts[t][a][b]*31.0*nSteps);
        if (a==50 && b==50) std::cout <<"Value is "<<avgStoreVals[t][a][b]<<std::endl;
      }
    }
  }


  //Get existing file info to copy to output file

  NcFile refFile(InputFiles[0].c_str());

  NcDim *inLatDim = refFile.get_dim("lat");
  NcDim *inLonDim = refFile.get_dim("lon");
  NcVar *inLatVar = refFile.get_var("lat");
  NcVar *inLonVar = refFile.get_var("lon");


//Output to file
  NcFile outfile(strfile_out.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);

  NcDim *outTime = outfile.add_dim("time", yearLen);
  DataVector<int> timeVals(yearLen);
  for (int t=0; t<yearLen; t++){
    timeVals[t] = t+1;
  }

  NcVar *outTimeVar = outfile.add_var("time", ncInt, outTime);
  outTimeVar->set_cur((long) 0);
  outTimeVar->put(&(timeVals[0]),yearLen);

  NcDim *outLat = outfile.add_dim("lat", nLat);
  NcDim *outLon = outfile.add_dim("lon", nLon);
  NcVar *outLatVar = outfile.add_var("lat", ncDouble, outLat);
  NcVar *outLonVar = outfile.add_var("lon", ncDouble, outLon);

  copy_dim_var(inLatVar,outLatVar);
  copy_dim_var(inLonVar,outLonVar);

  NcVar *outAvgIPV = outfile.add_var("AIPV", ncDouble, outTime, outLat, outLon);
  outAvgIPV->set_cur(0,0,0);
  outAvgIPV->put(&(avgStoreVals[0][0][0]),yearLen,nLat,nLon);

  refFile.close();
  outfile.close();
  }

  catch (Exception & e){
    std::cout<<e.ToString() <<std::endl;
  }
//  infile.close();
//  outfile.close();

  
}

