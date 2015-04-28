/////////////////////////////////////
///     \file moving_avg.cpp
///     \author Marielle Pinheiro
///     \version March 26, 2015

//#include "StitchBlobs.cpp"
#include "CLIVAR_block_utilities.h"
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

/*//Copied from StitchBlobs
void GetInputFileList(
        const std::string & strInputFileList,
                std::vector<std::string> & vecInputFiles
        ) {
                FILE * fp = fopen(strInputFileList.c_str(), "r");

                char szBuffer[1024];
                for (;;) {
                        fgets(szBuffer, 1024, fp);

                        if (feof(fp)) {
                                break;
                        }

                        // Remove end-of-line characters
                        for (;;) {
                                int nLen = strlen(szBuffer);
                                if ((szBuffer[nLen-1] == '\n') ||
                                        (szBuffer[nLen-1] == '\r') ||
                                        (szBuffer[nLen-1] == ' ')
                                ) {
                                        szBuffer[nLen-1] = '\0';
                                        continue;
                                }
                                break;
                        }

                        vecInputFiles.push_back(szBuffer);
        }

        if (vecInputFiles.size() == 0) {
                _EXCEPTION1("No files found in file \"%s\"", strInputFileList.c_str());
        }

        fclose(fp);
}
*/
//Copied from DetectCyclones
/*
void ParseTimeDouble(
        const std::string & strTimeUnits,
        const std::string & strTimeCalendar,
        double dTime,
        int & nDateYear,
        int & nDateMonth,
        int & nDateDay,
        int & nDateHour
) {
        // Get calendar type
        Time::CalendarType cal;
        if ((strTimeCalendar.length() >= 6) &&
                (strncmp(strTimeCalendar.c_str(), "noleap", 6) == 0)
        ) {
                cal = Time::CalendarNoLeap;

        } else if (
                (strTimeCalendar.length() >= 8) &&
                (strncmp(strTimeCalendar.c_str(), "standard", 8) == 0)
        ) {
                cal = Time::CalendarStandard;

        } else if (
                (strTimeCalendar.length() >= 9) &&
                (strncmp(strTimeCalendar.c_str(), "gregorian", 9) == 0)
        ) {
                cal = Time::CalendarStandard;

        } else {
                _EXCEPTION1("Unknown calendar type \"%s\"", strTimeCalendar.c_str());
        }

        // Time format is "days since ..."
        if ((strTimeUnits.length() >= 11) &&
            (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)
        ) {
                std::string strSubStr = strTimeUnits.substr(11);
                Time time(cal);
                time.FromFormattedString(strSubStr);

                int nDays = static_cast<int>(dTime);
                time.AddDays(nDays);

                int nSeconds = static_cast<int>(fmod(dTime, 1.0) * 86400.0);
                time.AddSeconds(nSeconds);

                Announce("Time (YMDS): %i %i %i %i",
                                time.GetYear(),
                                time.GetMonth(),
                                time.GetDay(),
                                time.GetSecond());

                nDateYear = time.GetYear();
                nDateMonth = time.GetMonth();
                nDateDay = time.GetDay();
                nDateHour = time.GetSecond() / 3600;

                //printf("%s\n", strSubStr.c_str());
        // Time format is "hours since ..."
        } else if (
            (strTimeUnits.length() >= 12) &&
            (strncmp(strTimeUnits.c_str(), "hours since ", 12) == 0)
        ) {
                std::string strSubStr = strTimeUnits.substr(12);
                Time time(cal);
                time.FromFormattedString(strSubStr);

                int nSeconds = static_cast<int>(fmod(dTime, 1.0) * 3600.0);
                time.AddSeconds(nSeconds);

                Announce("Time (YMDS): %i %i %i %i",
                                time.GetYear(),
                                time.GetMonth(),
                                time.GetDay(),
                                time.GetSecond());

                nDateYear = time.GetYear();
                nDateMonth = time.GetMonth();
                nDateDay = time.GetDay();
                nDateHour = time.GetSecond() / 3600;

                //printf("%s\n", strSubStr.c_str());

        } else {
                _EXCEPTIONT("Unknown \"time::units\" format");
        }
        //_EXCEPTION();
}
*/









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

  std::cout << "Opening first file."<<std::endl;

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

  double tRes = timeVec[1]-timeVec[0];

//  std::cout<<"Getting time unit information,"<<std::endl;

  NcAtt *attTime = timeVal->get_att("units");
  if (attTime == NULL){
    _EXCEPTIONT("Time variable has no units attribute.");
  }

  std::string strTimeUnits = attTime->as_string(0);

  NcAtt *attCal = timeVal->get_att("calendar");
  std::string strCalendar = attCal->as_string(0);

  std::cout<<"Time units: "<< strTimeUnits<<" Calendar: "<<strCalendar<<std::endl;


  
  int yearLen = 365;
  //if the calendar contains leap years, need to check for extra day
 // bool leap = false;
  int dateYear;
  int dateMonth;
  int dateDay;
  int dateHour;

  //Length of 31 days axis
  int arrLen = int(31.0/tRes);

  //3D matrix to store averaged values
  DataMatrix3D<double> avgStoreVals(yearLen,nLat,nLon);
  DataMatrix3D<double> avgCounts(yearLen,nLat,nLon);

  //Check whether file contains a Feb 29
  bool leap = false;

  int i=0;
  int leapYearIndex;

  //If there is a leap year, store the first occurrence of the leap year date
  while (i<nTime && leap==false){
    ParseTimeDouble(strTimeUnits, strCalendar, timeVec[i], dateYear,\
      dateMonth, dateDay, dateHour);
    if (dateMonth==2 && dateDay==29){
      std::cout<<"This file contains a leap year day at index "<<i<<std::endl;
      leap = true;
      leapYearIndex = i;
    }
    i++;
  }



  //Number of time steps per day
  int nSteps = 1/tRes;

  //Date index in which to start storing averaged values
  int dateIndex = int(timeVec[0])%yearLen+15;

  double endTime = timeVec[nTime-1];
  //Keeps track of index within 31-day array
  int currArrIndex = 0;

  //3D 31-day array
  DataMatrix3D<double> currFillData(arrLen,nLat,nLon);
//  std::cout<<"Initialized fill array."<<std::endl;

  bool newFile=false;
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
//  std::cout<<"Entering first while loop. tEnd is "<<tEnd<<\
    " and array length is "<<arrLen<<std::endl;
  //Initialize array with first 31 values
  while (currArrIndex<arrLen){
    for (int t=tStart; t<tEnd; t++){
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          if (leap){
            if (t>=leapYearIndex&&t<=(leapYearIndex+nSteps)){
              std::cout<<"leap year index is "<<leapYearIndex\
                <<" and t is currently "<<t<<". Skipping."<<std::endl;
              t++;
            }
          }
          currFillData[currArrIndex][a][b] = IPVData[t][a][b];
          std::cout<<"t is "<<t<<" and array index is "<<currArrIndex<<std::endl;
        }
      }
      currArrIndex+=1;
    }
    //if new file needs to be opened, open it
    if (currArrIndex < arrLen){
      infile.close();
      //Open new file
      x+=1;
      std::cout<<"Opening file "<<InputFiles[x]<<std::endl;
      NcFile infile(InputFiles[x].c_str());
      nTime = infile.get_dim("time")->size();
      nLat = infile.get_dim("lat")->size();
      nLon = infile.get_dim("lon")->size();

      //IPV variable
      NcVar *inPV = infile.get_var("IPV");
      DataMatrix3D<double> IPVData(nTime,nLat,nLon);

      inPV->set_cur(0,0,0);
      inPV->get(&(IPVData[0][0][0]),nTime,nLat,nLon);
      newFile = false;
      leap = false;

      //Time variable
      NcVar *timeVal = infile.get_var("time");
      DataVector<double> timeVec(nTime);
      timeVal->set_cur((long) 0);
      timeVal->get(&(timeVec[0]),nTime);

      std::cout<<"Check: end time was "<<endTime<<" and new start time is "\
        <<timeVec[0]<<std::endl;
      double contCheck = std::fabs(timeVec[0]-endTime);
      std::cout<<"Check: difference between old end time and new start time is "\
        <<contCheck<<" and tRes is "<< tRes<<std::endl;
      if (contCheck>tRes){
        _EXCEPTIONT("New file is not continuous with previous file."); 
      }

      //check for leap year days in this file
      i=0;
      while (i<nTime && leap==false){
        ParseTimeDouble(strTimeUnits, strCalendar, timeVec[i], dateYear,\
          dateMonth, dateDay, dateHour);
        if (dateMonth==2 && dateDay==29){
          std::cout<<"This file contains a leap year day at index "<<i<<std::endl;
          leap = true;
          leapYearIndex = i;
        }
        i++;
      }

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
  std::cout<<"Finished filling first array. Now averaging and saving to year array."<<std::endl;
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

  std::cout<<"Filled for date index "<<dateIndex<<std::endl;
  dateIndex+=1;
  currArrIndex = 0;

  //Enter while loop
  std::cout<<"Entering while loop!"<<std::endl;

  tStart = tEnd;
  tEnd = tStart + nSteps;
  //Check if new file needs to be opened
  if (tEnd >nTime){
    std::cout<<"Time dimension exceeded."<<std::endl;
    infile.close();
    x+=1;
    NcFile infile(InputFiles[x].c_str());
    std::cout<<"x is currently "<<x<<", opening file "<<InputFiles[x]<<std::endl;
    nTime = infile.get_dim("time")->size();
    nLat = infile.get_dim("lat")->size();
    nLon = infile.get_dim("lon")->size();

    //IPV variable
    NcVar *inPV = infile.get_var("IPV");
    DataMatrix3D<double> IPVData(nTime,nLat,nLon);

    inPV->set_cur(0,0,0);
    inPV->get(&(IPVData[0][0][0]),nTime,nLat,nLon);

    //Time variable
    NcVar *timeVal = infile.get_var("time");
    DataVector<double> timeVec(nTime);
    timeVal->set_cur((long) 0);
    timeVal->get(&(timeVec[0]),nTime);

    leap = false;
    i=0;
    while (i<nTime && leap==false){
      ParseTimeDouble(strTimeUnits, strCalendar, timeVec[i], dateYear,\
        dateMonth, dateDay, dateHour);
      if (dateMonth==2 && dateDay==29){
        std::cout<<"This file contains a leap year day at index "<<i<<std::endl;
        leap = true;
        leapYearIndex = i;
      }
      i++;
    }

    tStart = 0;
    tEnd = tStart + nSteps;
  }

  while (x<nFiles){
    //Advance to the next day and overwrite one day in the 31-day array
    std::cout<<"tStart: "<<tStart<<" tEnd: "<<tEnd<<std::endl;
    //Check that current replacement day isn't a Feb 29
    if (leap){
      if (tStart == leapYearIndex){
        tStart += nSteps;
        tEnd += nSteps;
        std::cout<<"Because of leap year, tStart/tEnd now "<<tStart<<" "<<tEnd<<std::endl;
      }
    }  


    std::cout<<"Replacing day starting at array index "<<currArrIndex<<std::endl;
    for (int t=tStart; t<tEnd; t++){
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          currFillData[currArrIndex][a][b] = IPVData[t][a][b];
        }
      }
      currArrIndex+=1;
    }
    std::cout<<"currArrIndex is currently "<<currArrIndex\
      <<" and array length is "<<arrLen<<std::endl;    
    //check for periodic boundary condition for 31 day array
    if (currArrIndex>=arrLen){
      std::cout<<"Periodic boundary: 31 day length exceeded."<<std::endl;
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
      }
    }
    std::cout<<"date index count is "<<avgCounts[dateIndex][1][1]<<std::endl;
    dateIndex+=1;
    //periodic boundary condition for year
    if(dateIndex>=yearLen){
      dateIndex-=yearLen;
      std::cout<<"dateIndex is now "<<dateIndex<<std::endl;
    }
    tStart = tEnd;
    tEnd = tStart + nSteps;
    std::cout<<"tEnd is currently "<<tEnd<<std::endl;

   //Check if new file needs to be opened
    if (tEnd>nTime){
      std::cout<<"Time dimension exceeded."<<std::endl;
      infile.close();
      std::cout<<"Closed file "<<InputFiles[x]<<std::endl;
      x+=1;
      if (x<nFiles){
      NcFile infile(InputFiles[x].c_str());
      std::cout<<"x is currently "<<x<<", opening file "<<InputFiles[x]<<std::endl;
      nTime = infile.get_dim("time")->size();
      nLat = infile.get_dim("lat")->size();
      nLon = infile.get_dim("lon")->size();

      //IPV variable
      NcVar *inPV = infile.get_var("IPV");
      DataMatrix3D<double> IPVData(nTime,nLat,nLon);

      inPV->set_cur(0,0,0);
      inPV->get(&(IPVData[0][0][0]),nTime,nLat,nLon);
      newFile = false;

      //Time variable
      NcVar *timeVal = infile.get_var("time");
      DataVector<double> timeVec(nTime);
      timeVal->set_cur((long) 0);
      timeVal->get(&(timeVec[0]),nTime);

      tStart = 0;
      tEnd = tStart + nSteps;}
    }

  }
  //average all values
  for (int t=0; t<yearLen; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        avgStoreVals[t][a][b] = avgStoreVals[t][a][b]/(avgCounts[t][a][b]*31.0);
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

