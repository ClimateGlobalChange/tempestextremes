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
      CommandLineString(varName, "varname","");
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
    if (varName == ""){
      _EXCEPTIONT("No variable name (--varname) specified");
    }

    //Multiplier for GH info
    double ghMult=1.;
    if (GHtoZ){
      ghMult = 1./9.8;
    }
    std::cout<<"ghmult is "<<ghMult<<std::endl;

    //STEP 1: GENERATE THE LIST OF FILES, GET UNITS AND SUCH
    //Create list of input files
    std::vector<std::string> InputFiles;
    GetInputFileList(fileList, InputFiles);
    int nFiles = InputFiles.size();

    //Open first file 
    NcFile reffile(InputFiles[0].c_str());
    std::cout<<"Reading in "<<InputFiles[0].c_str()<<" for reference"<<std::endl;
    int nTime = reffile.get_dim(tname.c_str())->size();
    int nLat = reffile.get_dim(latname.c_str())->size();
    int nLon = reffile.get_dim(lonname.c_str())->size();

    //time variable
    NcVar *timeVal = reffile.get_var(tname.c_str());
    if (timeVal==NULL){
      _EXCEPTIONT("Could not read time variable.");
    }
    DataVector<double> timeVec(nTime);
    timeVal->set_cur((long) 0);
    timeVal->get(&(timeVec[0]),nTime);

    //Time and Calendar attributes
    NcAtt *attCal = timeVal->get_att("calendar");
    std::string strCalendar;
    std::string strTimeUnits;
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
    double nt = 1./tRes;
    //Pressure axis (if 4D)
    int pIndex = 10000000;
    if (is4d){
      int nLev = reffile.get_dim(levname.c_str())->size();
      NcVar *levvar = reffile.get_var(levname.c_str());
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

    reffile.close();

    //Length of 61 days axis
    int arrLen = int(61.0/tRes);
    //Array time axis
    DataVector<double> arrTimeAxis(arrLen);
    DataVector<double> allTimes(arrLen);
    //Initialize with missing values
    for (int a=0; a<arrLen; a++){
        arrTimeAxis[a]=-999999.9;
    }
    //Fill Array
    DataMatrix3D<double> currFillData(arrLen,nLat,nLon);
    //Store matrix for summed values
    DataMatrix<double> sumValues(nLat,nLon);
    //Final year long array
    DataMatrix3D<double> smoothedValues(nDayYear,nLat,nLon);
    DataMatrix3D<double> smoothedCount(nDayYear,nLat,nLon);
    DataVector<double> fileTimeAxis(nDayYear);
    for (int t=0; t<nDayYear; t++){
        fileTimeAxis[t]=-999999.9;
    }
    //Keep track of array position
    int currArrIndex=0;
    //File days
    int currYear,currMonth,currDay,currHour;
    int nextYear,nextMonth,nextDay,nextHour;
    int prevYear,prevMonth,prevDay,prevHour;
    //Has missing values in any part of the array?
    bool hasMissing = false;
    int fillDay = 0;
    int yearDay = 0;
    double arrfrac = 1./(double(arrLen));
    //value of the time double for the middle of the 61 day array
    double midValue;
    //Keep track of how many days have passed so far
    int ntcount = 0;
    //Is this year a leap year?
    //bool isLeapYear = false;
    for (int x=0; x<nFiles; x++){
        NcFile infile(InputFiles[x].c_str());
        if (!infile.is_valid()){
            _EXCEPTION1("Could not open %s",InputFiles[x].c_str());
        }
        std::cout<<"Opened file "<<InputFiles[x].c_str()<<std::endl;
        NcVar * latvar = infile.get_var(latname.c_str());
        NcVar * lonvar = infile.get_var(lonname.c_str());
        NcVar * timeVal = infile.get_var(tname.c_str());
        if (timeVal==NULL){
            _EXCEPTIONT("Couldn't read time variable.");
        }   
        nTime = infile.get_dim(tname.c_str())->size();
        DataVector<double> timeVec(nTime);
        timeVal->set_cur(long(0));
        timeVal->get(&(timeVec[0]),nTime);
        //time units
        NcAtt *attTime = timeVal->get_att("units");
        if (attTime == NULL){
            _EXCEPTIONT("Time variable has no units attribute.");
        }
        strTimeUnits = attTime->as_string(0);
        //Read in the data
        NcVar * inData = infile.get_var(varName.c_str());
        if (inData==NULL){
            _EXCEPTIONT("Couldn't read in variable");
        }
        for (int t=0; t<nTime; t++){
            hasMissing=false;
            //Read in the time info
            ParseTimeDouble(strTimeUnits,strCalendar,timeVec[t],\
            currYear,currMonth,currDay,currHour);
            //The day that is being added to the averaging window
            fillDay = DayInYear(currMonth,currDay,strCalendar);
            //The day that is the center of the averaging window
            yearDay = fillDay - 31;
            //After first 30 days, dealing with periodic boundary
            if (ntcount>30*nt){
                if (yearDay<0){
                    yearDay+=365;
                }
            }
            //Fill the time axis with the number of days since the beginning of the year
            if (fileTimeAxis[yearDay]< -999999. && yearDay>=0){
                /*if (isLeapYear==true){
                  //From March 1st to 30th, need to subtract an extra day
                  if (currMonth==3 && currDay<31){
                    midValue-=tStep*nt;
                  }
                }*/
                //Fill with the time double minus 30 days
                //fileTimeAxis[yearDay] = midValue;
                fileTimeAxis[yearDay] = double(yearDay);
                std::cout<<"value "<<fileTimeAxis[yearDay]<<" added to "<<yearDay<<std::endl;
            }
            //The time double 
            midValue = timeVec[t]-(tStep*30*nt);
            //Special case when still filling in the first 30 days worth of data
            if (yearDay<0){
                midValue = -999999.9;
            }
            //std::cout<<"Adding fill for day "<<fillDay<<", which will go into day "<<yearDay<<std::endl;
            if (currMonth==2 && currMonth==29){
                //skip
                //isLeapYear = true;
            }else{
                //Get the data from the current time step in the netCDF
                if (is4d){
                    inData->set_cur(t,pIndex,0,0);
                    inData->get(&(currFillData[currArrIndex][0][0]),1,1,nLat,nLon);
                }else{
                    inData->set_cur(t,0,0);
                    inData->get(&(currFillData[currArrIndex][0][0]),1,nLat,nLon);
                }
                //Put the time value into the array time axis
                arrTimeAxis[currArrIndex]=timeVec[t];
                allTimes[currArrIndex]=timeVec[t];
                currArrIndex++;
                ntcount++;
                if (currArrIndex>=arrLen){
                    currArrIndex-=arrLen;
                }
                //Does the array time axis contain missing values?
                for (int z=0; z<arrLen; z++){
                    if (arrTimeAxis[z]< -999999.){
                        hasMissing=true;
                        break;
                    }
                }
                //If there are no missing values, sum up the data (makes a mean)
                sumValues.Initialize(nLat,nLon);
                if (hasMissing==false){
                    for (int z=0;z<arrLen;z++){
                        for (int a=0; a<nLat; a++){
                            for (int b=0; b<nLon; b++){
                                sumValues[a][b]+=currFillData[z][a][b]*arrfrac*ghMult;
                            }
                        }
                    }
                    //Now add that data to the year array (will be divided by the count at the end)
                    for (int a=0; a<nLat; a++){
                        for (int b=0; b<nLon; b++){
                            smoothedValues[yearDay][a][b]+=sumValues[a][b];
                            smoothedCount[yearDay][a][b]+=1.;
                        }
                    }

                }
            }
            //Check the current time step and the next time step to see if it's reached the end of the year
            ParseTimeDouble(strTimeUnits,strCalendar,midValue,\
            prevYear,prevMonth,prevDay,prevHour);
            ParseTimeDouble(strTimeUnits,strCalendar,midValue + tStep,\
            nextYear,nextMonth,nextDay,nextHour);
            std::cout<<"prevYear is "<<prevYear<<" and nextYear is "<<nextYear<<std::endl;
            //If it's reached the end of a year, write the file
            if (nextYear>prevYear){
                std::string outFileName = prefix + "_" + std::to_string(prevYear) + ".nc";
                std::string outTimeUnits = "days since " + std::to_string(prevYear);
                std::cout<<"inside while loop: writing file "<<outFileName<<std::endl;
                NcFile outfile(outFileName.c_str(),NcFile::Replace,NULL,0,NcFile::Offset64Bits);
                //Write the variables to the file
                NcDim * outTime = outfile.add_dim(tname.c_str(),nDayYear);
                NcDim * outlat = outfile.add_dim(latname.c_str(),nLat);
                NcDim * outlon = outfile.add_dim(lonname.c_str(),nLon);

                //Average the smoothed variable
                for (int z=0; z<nDayYear; z++){
                    for (int a=0; a<nLat; a++){
                        for (int b=0; b<nLon; b++){
                            smoothedValues[z][a][b] = smoothedValues[z][a][b]/smoothedCount[z][a][b];
                        }
                    }
                }
                //Add the time variable
                NcVar * outTimeVar = outfile.add_var(tname.c_str(),ncDouble,outTime);
                outTimeVar->set_cur(long(0));
                outTimeVar->put(&(fileTimeAxis[0]),nDayYear);
                std::string outCalendar;
                if (nDayYear<365){
                  outCalendar = "360_day";
                }else{
                  outCalendar = "noleap";
                }
                outTimeVar->add_att("units",outTimeUnits.c_str());
                outTimeVar->add_att("calendar",outCalendar.c_str());
                //lat/lon
                NcVar * outLatVar = outfile.add_var(latname.c_str(),ncDouble,outlat);
                NcVar * outLonVar = outfile.add_var(lonname.c_str(),ncDouble,outlon);
                copy_dim_var(latvar,outLatVar);
                copy_dim_var(lonvar,outLonVar);
                //Smoothed variable
                NcVar * outVar = outfile.add_var("SMOOTHED_Z500",ncDouble,outTime,outlat,outlon);
                outVar->set_cur(0,0,0);
                outVar->put(&(smoothedValues[0][0][0]),nDayYear,nLat,nLon);
                outfile.close();
                //Reset the time and smoothed variables
                for (int z=0; z<nDayYear;z++){
                    fileTimeAxis[z]=-999999.9;
                }
                std::cout<<"Reset file time axis"<<std::endl;
                smoothedValues.Initialize(nDayYear,nLat,nLon);
                smoothedCount.Initialize(nDayYear,nLat,nLon);
            }
        }
        infile.close();
    }
    //If it didn't reach the end of the year but there's more data to be written, write the last file
    //What was the last day in the year?
    if (yearDay<(nDayYear-1)){
      //open a reference file for the variables
      NcFile infile(InputFiles[nFiles-1].c_str());
      if (!infile.is_valid()){
        _EXCEPTION1("Could not open %s",InputFiles[nFiles-1].c_str());
      }
      //std::cout<<"Opened file "<<InputFiles[nFiles-1].c_str()<<std::endl;
      NcVar * latvar = infile.get_var(latname.c_str());
      NcVar * lonvar = infile.get_var(lonname.c_str());

      std::cout<<"Writing last file (did not fill the year out completely)"<<std::endl;
      std::cout<<"the year is "<<nextYear<<std::endl;
      std::string outFileName = prefix + "_" + std::to_string(prevYear) + ".nc";
      std::string outTimeUnits = "days since " + std::to_string(prevYear);

      NcFile outfile(outFileName.c_str(),NcFile::Replace,NULL,0,NcFile::Offset64Bits);
      //Write the variables to the file
      NcDim * outTime = outfile.add_dim(tname.c_str(),nDayYear);
      NcDim * outlat = outfile.add_dim(latname.c_str(),nLat);
      NcDim * outlon = outfile.add_dim(lonname.c_str(),nLon);

      for (int z=0; z<nDayYear; z++){
        if (fileTimeAxis[z]<-999999.){
          //fileTimeAxis[z]=fileTimeAxis[z-1]+(tStep*nt);
          fileTimeAxis[z] = fileTimeAxis[z-1] + 1.;
        }
      }
      //Average the smoothed variable
      for (int z=0; z<nDayYear; z++){
          for (int a=0; a<nLat; a++){
              for (int b=0; b<nLon; b++){
                  smoothedValues[z][a][b] = smoothedValues[z][a][b]/smoothedCount[z][a][b];
              }
          }
      }
      //Add the time variable
      NcVar * outTimeVar = outfile.add_var(tname.c_str(),ncDouble,outTime);
      outTimeVar->set_cur(long(0));
      outTimeVar->put(&(fileTimeAxis[0]),nDayYear);

      std::string outCalendar;
      if (nDayYear<365){
        outCalendar = "360_day";
      }else{
        outCalendar = "noleap";
      }
      outTimeVar->add_att("units",outTimeUnits.c_str());
      outTimeVar->add_att("calendar",outCalendar.c_str());
      //lat/lon
      NcVar * outLatVar = outfile.add_var(latname.c_str(),ncDouble,outlat);
      NcVar * outLonVar = outfile.add_var(lonname.c_str(),ncDouble,outlon);
      copy_dim_var(latvar,outLatVar);
      copy_dim_var(lonvar,outLonVar);
      //Smoothed variable
      NcVar * outVar = outfile.add_var("SMOOTHED_Z500",ncDouble,outTime,outlat,outlon);
      outVar->set_cur(0,0,0);
      outVar->put(&(smoothedValues[0][0][0]),nDayYear,nLat,nLon);
      outfile.close();
      infile.close();


    }


  }
  catch (Exception & e){
    std::cout<<e.ToString() <<std::endl;
  }
  
}