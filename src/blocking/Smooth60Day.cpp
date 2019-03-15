/////////////////////////////////////
///     \file Smooth61Day.cpp
///     \author Marielle Pinheiro
///     \version March 11, 2019


/*This code is modified from BlockingAvg.cpp and is
the intermediate step for detrending Z500 data. The
smoothing knocks down some of the higher frequency 
variability in the sub-daily Z500 data. The hope is 
to be able to use the CDO trend/subtrend tools in the
next step.  
*/
#include "BlockingUtilities.h"
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
    std::string varName;
    std::string avgName;
    std::string tname,latname,lonname,levname;
    bool missingFiles;
    std::string prefix;
    bool is4d;
    bool isHpa;
    bool GHtoZ;
    BeginCommandLine()
      CommandLineString(fileList, "inlist", "");
      CommandLineString(strfile_out, "out", "");
      CommandLineString(varName, "varname","");
      CommandLineString(avgName, "avgname","");
      CommandLineBool(missingFiles, "missing");
      CommandLineString(tname,"tname","time");
      CommandLineString(latname,"latname","lat");
      CommandLineString(lonname,"lonname","lon");
      CommandLineString(levname,"levname","lev");
      CommandLineString(prefix,"prefix","SMOOTHED");
      CommandLineBool(is4d,"is4D");
      CommandLineBool(isHpa,"hpa");
      CommandLineBool(GHtoZ,"g2z");
      ParseCommandLine(argc, argv);


    EndCommandLine(argv)
    AnnounceBanner();
    if (fileList == ""){
      _EXCEPTIONT("No file list (--inlist) provided");
    }
    if (strfile_out == ""){
      _EXCEPTIONT("No output file name (--out) provided");
    }
    if (varName == ""){
      _EXCEPTIONT("No variable name (--varname) specified");
    }
    if (avgName == ""){
      _EXCEPTIONT("No average name (--avgname) specified");
    }
    if (fileList == ""){
      _EXCEPTIONT("No file list (--inlist) specified");
    }

    //Multiplier for GH info
    double ghMult=1.;
    if (GHtoZ){
      ghMult = 1./9.8;
    }

    //Create list of input files
    std::vector<std::string> InputFiles;
    GetInputFileList(fileList, InputFiles);
    int nFiles = InputFiles.size();

    //Open first file 
    NcFile infile(InputFiles[0].c_str());
    std::cout<<"Reading in "<<InputFiles[0].c_str()<<std::endl;
    int nTime = infile.get_dim(tname.c_str())->size();
    int nLat = infile.get_dim(latname.c_str())->size();
    int nLon = infile.get_dim(lonname.c_str())->size();

    //time variable
    NcVar *timeVal = infile.get_var(tname.c_str());
    DataVector<double> timeVec(nTime);
    timeVal->set_cur((long) 0);
    timeVal->get(&(timeVec[0]),nTime);

    //Time and Calendar attributes
    NcAtt *attTime = timeVal->get_att("units");
    if (attTime == NULL){
      _EXCEPTIONT("Time variable has no units attribute.");
    }
    std::string strTimeUnits = attTime->as_string(0);

    NcAtt *attCal = timeVal->get_att("calendar");
    std::string strCalendar;
    if(attCal==NULL){
      std::cout<<"Calendar not found; setting to standard."<<std::endl;
      strCalendar = "standard";
    }else{
      strCalendar = attCal->as_string(0);
    }
    if (strncmp(strCalendar.c_str(), "gregorian",9)==0){
      strCalendar = "standard";
    }
    int nDayYear=365;
    if (strncmp(strCalendar.c_str(),"360_day",7)==0){
      nDayYear=360;
    }

    //Copy the lat and lon variables to data vectors for the output
    DataVector<double>latVec(nLat);
    DataVector<double>lonVec(nLon);
    NcVar * latvar = infile.get_var(latname.c_str());
    NcVar * lonvar = infile.get_var(lonname.c_str());

    //Resolution of the time axis
    double tRes;
    double tStep;
    tStep = timeVec[1]-timeVec[0];
    if ((strTimeUnits.length() >= 11) && \
      (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)){
      tRes = tStep;;
    }else if ((strTimeUnits.length()>= 14) &&\
      (strncmp(strTimeUnits.c_str(),"minutes since ",14)==0)){
      tRes = tStep/(24.*60.);
    }
    else {
      tRes = tStep/24.0;
    }
    //Length of 61 days axis
    int arrLen = int(61.0/tRes);

    //3D 61-day array
    DataMatrix3D<double> currFillData(arrLen,nLat,nLon);
    //IPV variable and data matrix
    NcVar *inPV = infile.get_var(varName.c_str());
    DataMatrix<double> IPVData(nLat, nLon);
    int pIndex = 10000000;
    if (is4d){
      int nLev = infile.get_dim(levname.c_str())->size();
      NcVar *levvar = infile.get_var(levname.c_str());
      DataVector<double> pVec(nLev);
      levvar->get(&pVec[0],nLev);

      double pval = 50000.0;
      if (isHpa){
        pval=500.;
      }
      for (int x=0; x<nLev; x++){
        if (std::fabs(pVec[x]-pval)<0.0001){
          pIndex = x;
	        break;
	      }
      }
      if (pIndex > 999999){
      	_EXCEPTIONT("Could not identify correct pressure level. Check file.");
      }
    }

    //File tracker
    int x=0;

    //File days
    int currYear = 0;
    int currMonth=0;
    int currDay=0;
    int currHour=0;

    //Begin adding data to the array. The first 61 days won't be smoothed!
    int currArrIndex=0;
    //Make a time axis for the store array
    DataVector<double>ArrTimeAxis(arrLen);
    DataVector<double>allTimes(arrLen);
    double eofTime=0;
    double nextTime=0;
    
    //1: Initially fill in data
    int currT = 0;
    std::string newFile = "FALSE";
    while (currArrIndex<arrLen){
      //NEW FILE IF BLOCK
      if (newFile=="TRUE"){
        currT = 0;
        //Open a new file
        NcFile infile(InputFiles[x].c_str());
        std::cout<<"Opening file "<<InputFiles[x].c_str()<<std::endl;
        if (!infile.is_valid()){
          _EXCEPTION1("Cannot open NetCDF file %s",InputFiles[x].c_str());
        }
        //get the new time dimension
        nTime = infile.get_dim(tname.c_str())->size();
        NcVar * timeVal = infile.get_var(tname.c_str());
        if (timeVal==NULL){
          _EXCEPTIONT("Couldn't read time variable.");
        }
        timeVec.Initialize(nTime);
        timeVal->set_cur(long(0));
        timeVal->get(&(timeVec[0]),nTime);

        NcVar * latvar = infile.get_var(latname.c_str());
        NcVar * lonvar = infile.get_var(lonname.c_str());
        //read in the new variable
        NcVar * inPV = infile.get_var(varName.c_str());
        if (inPV==NULL){
          _EXCEPTIONT("Couldn't read in variable");
        }
        newFile = "FALSE";
        if (x>0){
          //Is this file contiguous with the previous file?
          double nextTime = eofTime + tStep;
          ParseTimeDouble(strTimeUnits,strCalendar,\
          eofTime,currYear,currMonth,currDay,currHour);
          if (currMonth==2 && currDay==28 && nextTime<timeVec[0] ){
            int nsteps = 1/tRes;
            nextTime = eofTime + tStep*(nsteps+1);
          }
          if (timeVec[0]>nextTime){
            currArrIndex++;
            while(timeVec[0]>nextTime){
              ArrTimeAxis[currArrIndex]=-999999.9;
              allTimes[currArrIndex]=nextTime;
              nextTime+=tStep;
              currArrIndex++;
              std::cout<<"Filling time array axis with missing value at arr index="<<currArrIndex-1<<std::endl;
            }
          }
        }
      }
      if (!infile.is_valid()){
        _EXCEPTION1("Cannot read file %s",InputFiles[x].c_str());
      }
      //READ IN THE DATA
      ParseTimeDouble(strTimeUnits,strCalendar,timeVec[currT],\
        currYear,currMonth,currDay,currHour);
      if (currMonth==2 && currDay==29){
      //Skip leap days
        currT++;
      }
      else{
        if (is4d){
          inPV->set_cur(currT,pIndex,0,0);
          inPV->get(&(currFillData[currArrIndex][0][0]),1,1,nLat,nLon);
        }
        else{
          inPV->set_cur(currT,0,0);
          inPV->get(&(currFillData[currArrIndex][0][0]),1,nLat,nLon);
        }
        //Put the data into the fill index
        ArrTimeAxis[currArrIndex]=timeVec[currT];
        allTimes[currArrIndex]=timeVec[currT];
        currArrIndex++;
        currT++;
      }
      if (currT>=nTime){
        std::cout<<"About to close old file.";
        newFile="TRUE";
        eofTime = timeVec[nTime-1];
       // infile.close();
        x++;
      }
      
    }
    //2. Start the averaging
    //The start date is going to be at the start of the middle day
    int midDay = floor(arrLen/2);
    //Next while loop
    std::string hasMissing="FALSE";
    std::string writeFile="FALSE";
    double arrFrac = 1./double(arrLen);
    double arrTime;
    //reset the arr index 
    currArrIndex=0;
    //save the last time step and check if the year is going to switch in the next step
    double PrevDayDouble = allTimes[midDay];
    double NextDayDouble;
    DataMatrix3D<double> smoothedValues(nDayYear,nLat,nLon);
    DataMatrix3D<double> smoothedCount(nDayYear,nLat,nLon);
    DataVector<double> fileTimeAxis(nDayYear);
    for (int a=0; a<nDayYear; a++){
      //Fill the file time axis with missing values
      fileTimeAxis[a]=-999999.9;
    }

    //Next day values
    int NextYear, NextMonth, NextDay, NextHour;

    while(x<nFiles){
      //std::cout<<"current array index currently "<<currArrIndex<<" and currT is "<<currT<<", file "<<x+1<<" out of "<<nFiles<<std::endl;

      //Does a new file need to be opened?
      if (newFile=="TRUE"){
        currT = 0;
        //Open a new file
        NcFile infile(InputFiles[x].c_str());
        std::cout<<"Opening file "<<InputFiles[x].c_str()<<std::endl;
        if (!infile.is_valid()){
          _EXCEPTION1("Cannot open NetCDF file %s",InputFiles[x].c_str());
        }
        //get the new time dimension
        nTime = infile.get_dim(tname.c_str())->size();
        NcVar * timeVal = infile.get_var(tname.c_str());
        if (timeVal==NULL){
          _EXCEPTIONT("Couldn't read time variable.");
        }
        timeVec.Initialize(nTime);
        timeVal->set_cur(long(0));
        timeVal->get(&(timeVec[0]),nTime);
        //read in the new variable
        NcVar * inPV = infile.get_var(varName.c_str());
        if (inPV==NULL){
          _EXCEPTIONT("Couldn't read in variable");
        }
        newFile = "FALSE";
        if (x>0){
          //Is this file contiguous with the previous file?
          double nextTime = eofTime + tStep;
          ParseTimeDouble(strTimeUnits,strCalendar,\
            eofTime,currYear,currMonth,currDay,currHour);
          if (currMonth==2 && currDay==28 && nextTime<timeVec[0] ){
            double nsteps = 1./tRes;
            std::cout<<"There are "<<nsteps<<" timesteps per day"<<std::endl;
            nextTime = eofTime + tStep*(nsteps+1);
            std::cout<<"next time is now "<<nextTime<<std::endl;
          }
          if (timeVec[0]>nextTime){
            currArrIndex++;
            while(timeVec[0]>nextTime){
              allTimes[currArrIndex]=nextTime;
              ArrTimeAxis[currArrIndex]=-999999.9;
              nextTime+=tStep;
              currArrIndex++;
              NextDayDouble+=tStep;
              std::cout<<"Filling time array axis with missing value at arr index="<<currArrIndex-1<<std::endl;
            }
          }
        }
        std::cout<<"Finished reading in new file "<<InputFiles[x].c_str()<<std::endl;
      }
      //Initialize the matrix that will hope the smoothed value
      DataMatrix<double> sumValue(nLat,nLon);
      sumValue.Initialize(nLat,nLon);
      //std::cout<<"For t="<<currT<<" sumvalue is currently "<<sumValue[1][1]<<std::endl;
      //Do the sum of the 61-day array 
      //First check whether there are missing values
      for (int z=0; z<arrLen; z++){
        if (ArrTimeAxis[z]<-999999.){
          hasMissing = "TRUE";
          break;
        }
        else{
          hasMissing = "FALSE";
          //if there are no missing values, then add the 61-day 
          for (int a=0; a<nLat; a++){
            for (int b=0; b<nLon; b++){
              sumValue[a][b]+=currFillData[z][a][b]*arrFrac*ghMult;
            }
          }
        }
      }
      //std::cout<<"the array has been filled and sumvalue is now "<<sumValue[1][1]<<std::endl;
      //The data matrix has been filled
      //if there's a sum, copy that data matrix's sum to the year matrix
      if (hasMissing=="FALSE"){
        //Get the date from midDay
        //std::cout<<"mid day is "<<midDay<<std::endl;

        PrevDayDouble=allTimes[midDay];
        //std::cout<<"Prev day double is "<<PrevDayDouble<<" at index "<<midDay<<std::endl;
        NextDayDouble = PrevDayDouble + tStep;
        ParseTimeDouble(strTimeUnits,strCalendar,ArrTimeAxis[midDay],\
            currYear,currMonth,currDay,currHour);
        ParseTimeDouble(strTimeUnits,strCalendar,NextDayDouble,\
            NextYear,NextMonth,NextDay,NextHour);
        if (currYear != NextYear){
          writeFile = "TRUE";
        }else{
          writeFile = "FALSE";
        }
        //std::cout<<"Mid day date: "<<currYear<<"/"<<currMonth<<"/"<<currDay<<" at index "<<midDay<<std::endl;
        //Fill the appropriate day with the averaged value
        int dayIndex = DayInYear(currMonth,currDay,strCalendar)-1;
        //std::cout<<"day index is "<<dayIndex<<std::endl;
        //Fill in the time axis value
        if (fileTimeAxis[dayIndex]<-999999.){
          fileTimeAxis[dayIndex]=ArrTimeAxis[midDay];
        }else{
          //std::cout<<"Time axis value at index "<<dayIndex<<" is already set"<<std::endl;
        }

        //Add the values
        //std::cout<<"Smoothed values were previously "<<smoothedValues[dayIndex][1][1]<<" for day "<<dayIndex<<std::endl;
        for (int a=0; a<nLat; a++){
          for (int b=0; b<nLon; b++){
            smoothedValues[dayIndex][a][b] += sumValue[a][b];
            smoothedCount[dayIndex][a][b]+=1.;
          }
        }
        //std::cout<<"Smoothed values are now "<<smoothedValues[dayIndex][1][1]<<" at day index "<<dayIndex<<std::endl;
      }

      //Increment time step's value
      currT+=1;
     
      if (currT<nTime){
        //std::cout<<"CurrT is "<<currT<<std::endl;
        //READ IN THE DATA
        ParseTimeDouble(strTimeUnits,strCalendar,timeVec[currT],\
          currYear,currMonth,currDay,currHour);
        //std::cout<<"Reading in data for date "<<currYear<<currMonth<<currDay<<std::endl;
        if (currMonth==2 && currDay==29){
        //Skip leap days
          currT++;
        }
        else{
          if (is4d){
            inPV->set_cur(currT,pIndex,0,0);
            inPV->get(&(currFillData[currArrIndex][0][0]),1,1,nLat,nLon);
          }
          else{
            inPV->set_cur(currT,0,0);
            inPV->get(&(currFillData[currArrIndex][0][0]),1,nLat,nLon);
          }
        //Put the data into the fill index
          allTimes[currArrIndex]=timeVec[currT];
          ArrTimeAxis[currArrIndex]=timeVec[currT];
          //std::cout<<"Added value "<<timeVec[currT]<<"from "<<currT<<" to arr time axis at "<<currArrIndex<<" for time "<<currT<<std::endl;
          currArrIndex+=1;
          midDay+=1;
          if (midDay>=arrLen){
            midDay-=arrLen;
          }
          if(currArrIndex>=arrLen){
            currArrIndex-=arrLen;
          }
        }
      }
      if (currT>=nTime ){
        //std::cout<<"Curr t is "<<currT<<" and nTime is "<<nTime<<std::endl;
        std::cout<<"About to close old file.";

        newFile="TRUE";
        eofTime = timeVec[nTime-1];
        x++;
      }
      //Check the year of the next time step. If it's a different year, write the file
      if (writeFile=="TRUE"){
        //File name
        std::string outFileName = prefix + "_" + std::to_string(currYear-1) + ".nc";
        std::cout<<"inside while loop: writing file "<<outFileName<<std::endl;
        NcFile outfile(outFileName.c_str(),NcFile::Replace,NULL,0,NcFile::Offset64Bits);
        //Write the variables to the file
        NcDim * outTime = outfile.add_dim(tname.c_str(),nDayYear);
        NcDim * outlat = outfile.add_dim(latname.c_str(),nLat);
        NcDim * outlon = outfile.add_dim(lonname.c_str(),nLon);

        for (int t=0; t<nDayYear; t++){
          for (int a=0; a<nLat; a++){
            for (int b=0; b<nLon; b++){
              smoothedValues[t][a][b]=smoothedValues[t][a][b]/smoothedCount[t][a][b];
            }
          }
        }
        //Add the time variable
        //Scan for the first instance of non-missing time
        //if there are missing values, work back and fill them in
        int breakIndex=0;
        for (int a=0; a<nDayYear; a++){
          if (fileTimeAxis[a]< -999999.){
            //Missing, move on
          }else{
            std::cout<<"Found break point at index "<<a<<std::endl;
            breakIndex=a;
            break;
          }
        }
        if (breakIndex>0){
          //Work backwards and fill in the double values
          for (int a=breakIndex; a>=0; a--){
            if (fileTimeAxis[a]<-999999.){
              fileTimeAxis[a]=fileTimeAxis[a+1]-(tStep/tRes);
              //std::cout<<"Replaced missing value with "<<fileTimeAxis[a]<<" at "<<a<<std::endl;
            }
          }
        }
        //Now check forward as well
        for (int a=breakIndex; a<nDayYear; a++){
          if (fileTimeAxis[a]< -999999.){
            fileTimeAxis[a]=fileTimeAxis[a-1]+(tStep/tRes);
            //std::cout<<"going forward, replaced missing value with "<<fileTimeAxis[a]<<" at "<<a<<std::endl;
          }
        }
        //Write the time axis to file
        NcVar * outTimeVar = outfile.add_var(tname.c_str(),ncDouble,outTime);
        outTimeVar->set_cur(long(0));
        outTimeVar->put(&(fileTimeAxis[0]),nDayYear);
        outTimeVar->add_att("units",strTimeUnits.c_str());
        outTimeVar->add_att("calendar",strCalendar.c_str());
        //Write lat/lon axes
        NcVar * outLatVar = outfile.add_var(latname.c_str(),ncDouble,outlat);
        NcVar * outLonVar = outfile.add_var(lonname.c_str(),ncDouble,outlon);
        copy_dim_var(latvar,outLatVar);
        copy_dim_var(lonvar,outLonVar);
        //Write the smoothed variable
        NcVar * outSmoothedVar = outfile.add_var("SMOOTHED_Z500",ncDouble,outTime,outlat,outlon);
        outSmoothedVar->set_cur(0,0,0);
        outSmoothedVar->put(&(smoothedValues[0][0][0]),nDayYear,nLat,nLon);


        //reset the time axis values for the next go round
        for (int a=0; a<nDayYear; a++){
          //Fill the file time axis with missing values
          fileTimeAxis[a]=-999999.9;
        }
        //smoothedValues.Initialize(nDayYear,nLat,nLon);
        std::cout<<"Resetting smoothed values"<<std::endl;
        for (int t=0; t<nDayYear; t++){
          for (int a=0; a<nLat; a++){
            for (int b=0; b<nLon; b++){
              smoothedValues[t][a][b]=0.;
              smoothedCount[t][a][b]=0.;
            }
          }
        }
      }
    }
    //Write the last file for the remaining data
    //File name
    std::string outFileName = prefix + "_" + std::to_string(currYear) + ".nc";
    std::cout<<"outside while loop: writing file "<<outFileName<<std::endl;
    NcFile outfile(outFileName.c_str(),NcFile::Replace,NULL,0,NcFile::Offset64Bits);
    //Write the variables to the file
    NcDim * outTime = outfile.add_dim(tname.c_str(),nDayYear);
    NcDim * outlat = outfile.add_dim(latname.c_str(),nLat);
    NcDim * outlon = outfile.add_dim(lonname.c_str(),nLon);

    for (int t=0; t<nDayYear; t++){
      for (int a=0; a<nLat; a++){
        for (int b=0; b<nLon; b++){
          smoothedValues[t][a][b]=smoothedValues[t][a][b]/smoothedCount[t][a][b];
        }
      }
    }
    //Add the time variable
    //Scan for the first instance of non-missing time
    //if there are missing values, work back and fill them in
    int breakIndex=0;
    for (int a=0; a<nDayYear; a++){
      //std::cout<<"file axis is "<<fileTimeAxis[a]<<" at index "<<a<<std::endl;
      if (fileTimeAxis[a]< -999999.){
        //Missing, move on
      }else{
        //std::cout<<"Found break point at index "<<a<<std::endl;
        breakIndex=a;
        break;
      }
    }
    if (breakIndex>0){
      //Work backwards and fill in the double values
      for (int a=breakIndex; a>=0; a--){
        if (fileTimeAxis[a]<-999999.){
          fileTimeAxis[a]=fileTimeAxis[a+1]-(tStep/tRes);
          //std::cout<<"Replaced missing value with "<<fileTimeAxis[a]<<" at "<<a<<std::endl;
        }
      }
    }
    //Now check forward as well
    for (int a=breakIndex; a<nDayYear; a++){
      if (fileTimeAxis[a]< -999999.){
        fileTimeAxis[a]=fileTimeAxis[a-1]+(tStep/tRes);
        //std::cout<<"going forward, replaced missing value with "<<fileTimeAxis[a]<<" at "<<a<<std::endl;
      }
    }
    //Write the time axis to file
    NcVar * outTimeVar = outfile.add_var(tname.c_str(),ncDouble,outTime);
    outTimeVar->set_cur(long(0));
    outTimeVar->put(&(fileTimeAxis[0]),nDayYear);
    outTimeVar->add_att("units",strTimeUnits.c_str());
    outTimeVar->add_att("calendar",strCalendar.c_str());
    //Write lat/lon axes
    NcVar * outLatVar = outfile.add_var(latname.c_str(),ncDouble,outlat);
    NcVar * outLonVar = outfile.add_var(lonname.c_str(),ncDouble,outlon);
    copy_dim_var(latvar,outLatVar);
    copy_dim_var(lonvar,outLonVar);
    //Write the smoothed variable
    NcVar * outSmoothedVar = outfile.add_var("SMOOTHED_Z500",ncDouble,outTime,outlat,outlon);
    outSmoothedVar->set_cur(0,0,0);
    outSmoothedVar->put(&(smoothedValues[0][0][0]),nDayYear,nLat,nLon);
    //std::cout<<"At time t=90, the value is "<<smoothedValues[90][1][1]<<std::endl;


  }
  catch (Exception & e){
    std::cout<<e.ToString() <<std::endl;
  }
  
}

