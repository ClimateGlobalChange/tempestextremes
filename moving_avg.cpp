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

//Copied from StitchBlobs
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

/*bool containsLeap(NcVar *timeVals, std::string units){
  std::vector<std::string> unitsSplit;
  std::stringstream unitsStream(units);
  std::string word;
  while (std::getline(unitsStream, word, " ")){
    unitsSplit.push_back(word);
  }
  std::cout<<"Test: first word is "<<unitsSplit[0]<<std::endl;
}*/

int main(int argc, char **argv){
  NcError error(NcError::verbose_nonfatal);

  try{
  std::string fileList;
  bool leap;
  std::string strfile_out;

  BeginCommandLine()
    CommandLineString(fileList, "in", "");
    CommandLineString(strfile_out, "out", "");
    CommandLineBool(leap, "leap");
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

  double tRes = timeVec[1]-timeVec[0];

  int yearLen = 365;
  //if the calendar contains leap years, need to check for extra day
  if (leap){
    NcAtt *timeUnits = timeVal->get_att("units");
    if (timeUnits==NULL){
      _EXCEPTIONT("Time variable has no units attribute.");
    }
    std::string strUnits = timeUnits->as_string(0);
   // NcAtt *timeCal = timeVal->get_att("calendar");
    std::cout<<"Time units are "<< strUnits <<std::endl;
   
  }


  //Length of 31 days axis
  int arrLen = int(31.0/tRes);

  //3D matrix to store averaged values
  DataMatrix3D<double> avgStoreVals(yearLen,nLat,nLon);
  DataMatrix3D<double> avgCounts(yearLen,nLat,nLon);

  //Date in which to start storing averaged values
  int dateIndex = int(timeVec[0])%yearLen+15;
  //Start date of first file
  double startTime = std::fmod(timeVec[0],yearLen);
  std::cout<<"Will start storing value "<<startTime<<" at index "<<dateIndex<<std::endl;

  //Keeps track of index within 31-day array
  int currArrIndex = 0;

  //number of time steps in day
  int nSteps = 1/tRes;
  std::cout<<"Number of steps in day is "<<nSteps<<std::endl;

  //3D 31-day array
  DataMatrix3D<double> currFillData(arrLen,nLat,nLon);
  std::cout<<"Initialized fill array."<<std::endl;

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
  std::cout<<"Entering first while loop. tEnd is "<<tEnd<<\
    " and array length is "<<arrLen<<std::endl;
  //Initialize array with first 31 values
  while (currArrIndex<arrLen){
    for (int t=tStart; t<tEnd; t++){
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          currFillData[currArrIndex][a][b] = IPVData[t][a][b];
        }
      }
      currArrIndex+=1;
      //currTimeIndex+=1;
    }
   // std::cout<<"Current fill array index is "<<currArrIndex<<\
      " and time index is "<<currTimeIndex<<std::endl;
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

      //Time variable
      NcVar *timeVal = infile.get_var("time");
      DataVector<double> timeVec(nTime);
      timeVal->set_cur((long) 0);
      timeVal->get(&(timeVec[0]),nTime);
      //currTimeIndex = 0;

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
      //avgStoreVals[dateIndex][a][b] = avgStoreVals[dateIndex][a][b]/arrLen;
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

    tStart = 0;
    tEnd = tStart + nSteps;
  }

  while (x<nFiles){
    //Advance to the next day and overwrite one day in the 31-day array
    std::cout<<"tStart: "<<tStart<<" tEnd: "<<tEnd<<std::endl;

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

