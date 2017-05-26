///////////////////////////////////////////////////
///
///           \file BlockingZThresh.cpp
///           \author Marielle Pinheiro
///           \date May 8, 2017
///
//////////////////////////////////////////////////

/*
This file calculates the threshold value at each (x,y)
for the corresponding day in the year. The threshold is
defined as 1.5 times the standard deviation of Z500 values.

It takes an input list of instantaneous Z500 values, as well
as the averaged file, and calculates the standard deviation 
of the instantaneous values. The output threshold is saved to
the average file.
*/
#include "netcdfcpp.h"
#include "DataVector.h"
#include "DataMatrix3D.h"
#include "DataMatrix.h"
#include "Announce.h"
#include "CommandLine.h"
#include "Exception.h"
#include "BlockingUtilities.h"
#include "NetCDFUtilities.h"
#include "DFT.h"

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <complex>

int main(int argc, char ** argv){
  try{
    std::string outFile,fileList,avgName,avgFile,varName,tname,latname,lonname;
    BeginCommandLine()
      CommandLineString(outFile,"outfile","");
      CommandLineString(fileList,"inlist","");
      CommandLineString(avgFile,"avgfile","");
      CommandLineString(avgName,"avgname","");
      CommandLineString(varName,"varname","");
      CommandLineString(tname,"tname","time");
      CommandLineString(latname,"latname","lat");
      CommandLineString(lonname,"lonname","lon");
      ParseCommandLine(argc,argv);
    EndCommandLine(argv)
    AnnounceBanner();

    if (fileList == ""){
      _EXCEPTIONT("No file list provided.");
    }
		if (avgFile==""){
			_EXCEPTIONT("No average file name provided.");
		}
    if (varName==""){
      _EXCEPTIONT("Need to provide input variable name (--varname).");
    }
    if (avgName == ""){
      _EXCEPTIONT("Need to provide input average variable name (--avgname).");
    }
    if (outFile==""){
      _EXCEPTIONT("Need to provide output file name.");
    }
    std::vector<std::string> vecFiles;
    GetInputFileList(fileList,vecFiles);
    int nFiles = vecFiles.size();
		
    //Open up the averages file
    NcFile avgin(avgFile.c_str());
    if (!avgin.is_valid()){
      _EXCEPTION1("Unable to open file %s for reading",avgFile.c_str());
    }
    //axis lengths
    int latLen,lonLen,dLen;
  
    NcDim *latDim = avgin.get_dim(latname.c_str());
    latLen = latDim->size();
    NcDim *lonDim = avgin.get_dim(lonname.c_str());
    lonLen = lonDim->size();
    NcDim *dDim = avgin.get_dim(tname.c_str());
    dLen = dDim->size();
		
    //Read in the matrix of average Z500 values
    NcVar *avgVar = avgin.get_var(avgName.c_str());
    DataMatrix3D<double> avgMat(dLen,latLen,lonLen);
    avgVar->set_cur(0,0,0);
    avgVar->get(&(avgMat[0][0][0]),dLen,latLen,lonLen);
				
    //Open up the first file
    NcFile readin(vecFiles[0].c_str());
    if (!readin.is_valid()){
      _EXCEPTION1("Unable to open file %s for reading",vecFiles[0].c_str());
    }
    int tLen;
    tLen = readin.get_dim(tname.c_str())->size();
    //Initialize storage array: day, lat, lon
    DataMatrix3D<double> storeMat(dLen,latLen,lonLen);
    //Initialize counts array: day, lat, lon
    DataMatrix3D<double> countsMat(dLen, latLen, lonLen);
    //Initialize input data: time, lat, lon
    DataMatrix3D<double> inputData(tLen,latLen,lonLen);
    DataVector<double> timeVec(tLen);

    //Close file (then re-open in loop)
    readin.close();
		
    //Loop through all available files and store values in appropriate day
    for (int x=0; x<nFiles; x++){
      NcFile readin(vecFiles[x].c_str());
      if (!readin.is_valid()){
        _EXCEPTION1("Unable to open file %s for reading",vecFiles[x].c_str());
      }

      tLen = readin.get_dim(tname.c_str()) ->size();
      //std::cout<<"Reading in file #"<<x<<", "<<vecFiles[x].c_str()<<std::endl;
      //Input data
      inputData.Initialize(tLen, latLen, lonLen);
      NcVar *inputVar = readin.get_var(varName.c_str());
      inputVar->set_cur(0,0,0);
      inputVar->get(&(inputData[0][0][0]),tLen,latLen,lonLen);

      //Time and calendar
      timeVec.Initialize(tLen);
      NcVar *timeVar = readin.get_var(tname.c_str());
      timeVar->set_cur((long)0);
      timeVar->get(&(timeVec[0]),tLen);

      NcAtt *attTime = timeVar->get_att("units");
      if (attTime==NULL){
        _EXCEPTIONT("Time variable has no units attribute.");
      }
      std::string strTimeUnits = attTime->as_string(0);

      NcAtt *attCal = timeVar->get_att("calendar");
      if (attCal==NULL){
        _EXCEPTIONT("Time variable has no calendar attribute.");
      }
      std::string strCalendar = attCal->as_string(0);
      if (strncmp(strCalendar.c_str(),"gregorian",9)==0){
        strCalendar = "standard";
      }

      int dateYear, dateMonth, dateDay, dateHour, dayIndex;
      bool leap;
			double avgValue,currDev;
      //std::cout<<"Storing values."<<std::endl;
      //Store the values and counts in the matrices at the appropriate day index
      for (int t=0; t<tLen; t++){
        ParseTimeDouble(strTimeUnits, strCalendar, timeVec[t],\
          dateYear, dateMonth, dateDay, dateHour);
        leap = checkFileLeap(strTimeUnits, strCalendar, dateYear,\
          dateMonth, dateDay, dateHour, timeVec[t]);
        if (leap == false){
          dayIndex = DayInYear(dateMonth, dateDay)-1;
          //std::cout<<"Day index is "<<dayIndex<<std::endl;
          for (int a=0; a<latLen; a++){
            for (int b=0; b<lonLen; b++){
              avgValue = avgMat[dayIndex][a][b];
              currDev = inputData[t][a][b]-avgValue;
              storeMat[dayIndex][a][b]+= (currDev*currDev);
              countsMat[dayIndex][a][b]+= 1.;
            }
          }
        }
      }
      readin.close();
    }
    //Calculating 1.5*SD
    for (int d=0;d<dLen;d++){
      for (int a=0;a<latLen;a++){
        for (int b=0;b<lonLen;b++){
          storeMat[d][a][b]=(1.5*std::sqrt(storeMat[d][a][b]/countsMat[d][a][b]));
	}
      }
    }
		
    //Check if the NetCDF variable already exists in the file
 
    /*std::cout<<"About to write new variables"<<std::endl;
    NcBool exists_var = avgin.get_var("THRESHOLD_Z")->is_valid();
    std::cout<<"The status of the variable is "<<exists_var<<std::endl;
    NcVar *checkThresh;
    if (avgin.get_var("THRESHOLD_Z")->is_valid()==false){
      std::cout<<"Adding new threshold variable."<<std::endl;
      checkThresh = avgin.add_var("THRESHOLD_Z",ncDouble,dDim,latDim,lonDim);      
    }else{
      checkThresh = avgin.get_var("THRESHOLD_Z");
    }*/

    NcFile fileout(outFile.c_str(),NcFile::Replace,NULL,0,NcFile::Offset64Bits);
    NcDim *outTime = fileout.add_dim(tname.c_str(),dLen);
    NcDim *outLat = fileout.add_dim(latname.c_str(),latLen);
    NcDim *outLon = fileout.add_dim(lonname.c_str(),lonLen);
    NcVar *outTvar = fileout.add_var(tname.c_str(),ncDouble,outTime);
    NcVar *outLatvar = fileout.add_var(latname.c_str(),ncDouble,outLat);
    NcVar *outLonvar = fileout.add_var(lonname.c_str(),ncDouble,outLon);

    copy_dim_var(avgin.get_var(tname.c_str()),outTvar);
    copy_dim_var(avgin.get_var(latname.c_str()),outLatvar);
    copy_dim_var(avgin.get_var(lonname.c_str()),outLonvar);

    NcVar *checkThresh = fileout.add_var("THRESHOLD",ncDouble,outTime,outLat,outLon);
    checkThresh->set_cur(0,0,0);
    checkThresh->put(&(storeMat[0][0][0]),dLen,latLen,lonLen);		

    //Just for fun... the DFT of the SDs
    std::vector<double> inputDaily(dLen);
    std::vector<std::complex<double> >FourierCoefs(dLen);
    std::vector<double>outputDaily(dLen);
    DataMatrix3D<double> outputMat(dLen,latLen,lonLen);

    for (int a=0; a<latLen; a++){
      for (int b=0; b<lonLen; b++){
        for (int d=0; d<dLen; d++){
          inputDaily[d] = storeMat[d][a][b];
        }
        FourierCoefs = DFT(inputDaily,6);
        outputDaily = IDFT(FourierCoefs);
        for (int d=0; d<dLen; d++){
          outputMat[d][a][b] = outputDaily[d];
        }
      }
    }
    
    NcVar *checkDFTThresh = fileout.add_var("THRESHOLD_DFT",ncDouble,outTime,outLat,outLon);
  /*  if (avgin.get_var("THRESHOLD_Z_DFT")->is_valid()==false){
      std::cout<<"Adding new DFT threshold variable."<<std::endl;
      checkDFTThresh = avgin.add_var("THRESHOLD_Z_DFT",ncDouble,dDim,latDim,lonDim);
    }
    else{
      checkDFTThresh = avgin.get_var("THRESHOLD_Z_DFT");
    }*/
    checkDFTThresh->set_cur(0,0,0);
    checkDFTThresh->put(&(outputMat[0][0][0]),dLen,latLen,lonLen);

    avgin.close();
    fileout.close();
  }
  catch (Exception &e){
    std::cout<<e.ToString()<<std::endl;
  }
}
