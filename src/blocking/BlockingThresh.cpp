///////////////////////////////////////////////////
///
///           \file BlockingThresh.cpp
///           \author Marielle Pinheiro
///           \date May 8, 2017
///
//////////////////////////////////////////////////

/*
This file calculates the threshold value at each (x,y)
for the corresponding day in the year. The threshold is
defined as 1.5 times the standard deviation of Z500 values.

It takes an input list of instantaneous deviation  values, as well
as the average of the deviations file, and calculates the standard deviation 
of the instantaneous deviation values.
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
  NcError error(NcError::silent_nonfatal);
  try{
    std::string outFile,fileList,avgName,avgFile,varName,tname,latname,lonname;
    int nWavesLat,nWavesLon,nWavesTime;
    BeginCommandLine()
      CommandLineString(outFile,"outfile","");
      CommandLineString(fileList,"inlist","");
      CommandLineString(avgFile,"davgfile","");
      CommandLineString(avgName,"davgname","");
      CommandLineString(varName,"dvarname","");
      CommandLineString(tname,"tname","time");
      CommandLineString(latname,"latname","lat");
      CommandLineString(lonname,"lonname","lon");
      CommandLineInt(nWavesLat,"nLat",2);
      CommandLineInt(nWavesLon,"nLon",2);
      CommandLineInt(nWavesTime,"nTime",8);
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

    //Open lats vector
    NcVar *latRef = avgin.get_var(latname.c_str());
    DataVector <double> latVec(latLen);
    latRef->set_cur((long) 0);
    latRef->get(&(latVec[0]),latLen);
		
    //Read in the matrix of average values
    NcVar *avgVar = avgin.get_var(avgName.c_str());
    DataMatrix3D<double> avgMat(dLen,latLen,lonLen);
    avgVar->set_cur(0,0,0);
    avgVar->get(&(avgMat[0][0][0]),dLen,latLen,lonLen);
    std::cout<<"Opening input files"<<std::endl;	

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
    DataMatrix<double> inputData(latLen,lonLen);
    DataVector<double> timeVec(tLen);

    //Close file (then re-open in loop)
    readin.close();
		
    //Loop through all available files and store values in appropriate day
    //Part 1: Calculating the average
    for (int x=0; x<nFiles; x++){
      NcFile readin(vecFiles[x].c_str());
      if (!readin.is_valid()){
        _EXCEPTION1("Unable to open file %s for reading",vecFiles[x].c_str());
      }

      tLen = readin.get_dim(tname.c_str()) ->size();
      //std::cout<<"Reading in file #"<<x<<", "<<vecFiles[x].c_str()<<std::endl;
      //Input data
    //  inputData.Initialize(tLen, latLen, lonLen);
      NcVar *inputVar = readin.get_var(varName.c_str());
    //  inputVar->set_cur(0,0,0);
    //  inputVar->get(&(inputData[0][0][0]),tLen,latLen,lonLen);

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
      std::string strCalendar;
      NcAtt *attCal = timeVar->get_att("calendar");
      if (attCal==NULL){
        strCalendar = "standard";
      }else{
        strCalendar = attCal->as_string(0);
      }
      if (strncmp(strCalendar.c_str(),"gregorian",9)==0){
        strCalendar = "standard";
      }

      int dateYear, dateMonth, dateDay, dateHour, dayIndex;
      bool leap;
			double avgValue,currDev;
      //std::cout<<"Storing values."<<std::endl;
      //Store the values and counts in the matrices at the appropriate day index
      for (int t=0; t<tLen; t++){
        inputData.Initialize(latLen,lonLen);
        inputVar->set_cur(t,0,0);
        inputVar->get(&(inputData[0][0]),1,latLen,lonLen);
        ParseTimeDouble(strTimeUnits, strCalendar, timeVec[t],\
          dateYear, dateMonth, dateDay, dateHour);
        leap = checkFileLeap(strTimeUnits, strCalendar, dateYear,\
          dateMonth, dateDay, dateHour, timeVec[t]);
        if (leap == false){
          dayIndex = DayInYear(dateMonth, dateDay)-1;
          for (int a=0; a<latLen; a++){
            for (int b=0; b<lonLen; b++){
              avgValue = avgMat[dayIndex][a][b];
              currDev = inputData[a][b]-avgValue;
              storeMat[dayIndex][a][b]+= (currDev*currDev);
              countsMat[dayIndex][a][b]+= 1.;
           /*   if (a==70 && b==300){
              std::cout<<"DEBUG: t is "<<t<<" and dayIndex is "<<dayIndex<<", avg value is "\
                <<avgValue<<", dev is "<<currDev<<", just added "<<currDev*currDev<<\
                ", countsMat is now "<<countsMat[dayIndex][a][b]<<std::endl;
              }*/
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
        /*  if (a==70 && b==300){
          std::cout<<"DEBUG: for d="<<d<<", storeMat is "<<storeMat[d][a][b]<<std::endl;
          }*/
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
/*
    NcVar *counts = fileout.add_var("COUNTS",ncDouble,outTime,outLat,outLon);
    counts->set_cur(0,0,0);
    counts->put(&(countsMat[0][0][0]),dLen,latLen,lonLen);
*/
 /*   NcVar *checkThresh = fileout.add_var("THRESHOLD",ncDouble,outTime,outLat,outLon);
    checkThresh->set_cur(0,0,0);
    checkThresh->put(&(storeMat[0][0][0]),dLen,latLen,lonLen);		
*/

    std::cout<<"calculating DFT of threshold"<<std::endl;
    //First, try zonal average
    std::vector<double>zonalDaily(lonLen);
    std::vector<std::complex<double> >FC(lonLen);
    std::vector<double> zonalOut(lonLen);
    DataMatrix3D<double> zonalMat(dLen,latLen,lonLen);

    for (int d=0; d<dLen; d++){
      for (int a=0; a<latLen; a++){
        for (int b=0; b<lonLen; b++){
          zonalDaily[b] = storeMat[d][a][b];
        }
        FC = DFT(zonalDaily,nWavesLon);
        zonalOut = IDFT(FC);
        for (int b=0; b<lonLen; b++){
          zonalMat[d][a][b] = zonalOut[b];
        }
      }
    }
    //Add on a lat average
    //Do separate averages for NH and SH
    //Find index for equator
    int eqIndex;
    for (int e=0; e<latLen; e++){
      if (std::fabs(latVec[e])<0.0001){
        eqIndex = e;
        break;
      }
    }
    int lat1 = latVec[0];
    //Case 1: lat vector goes from NH to SH
    int NHLen,SHLen;
    if (lat1>0){
      NHLen = eqIndex;
      SHLen = latLen - eqIndex;
    }
    else{
      NHLen = latLen - eqIndex;
      SHLen = eqIndex;
    }

    std::vector<double>NH(NHLen);
    std::vector<double>SH(SHLen);
    std::vector<std::complex<double> > FCNH(NHLen);
    std::vector<std::complex<double> > FCSH(SHLen);
    std::vector<double> NHout(NHLen);
    std::vector<double> SHout(SHLen);
    int a1;
    DataMatrix3D<double>zmMat(dLen,latLen,lonLen);
    for (int d=0; d<dLen; d++){
      for (int b=0; b<lonLen; b++){
        for (int a=0; a<eqIndex; a++){
          NH[a] = zonalMat[d][a][b];
        }
        for (int a=eqIndex; a<latLen; a++){
          a1 = a-eqIndex;
          SH[a1] = zonalMat[d][a][b];
        }
        FCNH = DFT(NH,nWavesLat);
        NHout = IDFT(FCNH);
        FCSH = DFT(SH,nWavesLat);
        SHout = IDFT(FCSH);
        for (int a=0; a<eqIndex; a++){
          zmMat[d][a][b] = NHout[a];
        }
        for (int a=eqIndex; a<latLen; a++){
          a1 = a-eqIndex;
          zmMat[d][a][b] = SHout[a1];
        }
      }
    }

    //Alternate attempt: moving average (if DFT not appropriate)


/*
    NcVar * zmThresh = fileout.add_var("THRESHOLD_AVG",ncDouble,outTime,outLat,outLon);
    zmThresh->set_cur(0,0,0);
    zmThresh->put(&(zmMat[0][0][0]),dLen,latLen,lonLen);
  */

  //Just for fun... the DFT of the SDs
    std::vector<double> inputDaily(dLen);
    std::vector<std::complex<double> >FourierCoefs(dLen);
    std::vector<double>outputDaily(dLen);
    DataMatrix3D<double> outputMat(dLen,latLen,lonLen);


    for (int a=0; a<latLen; a++){
      for (int b=0; b<lonLen; b++){
        for (int d=0; d<dLen; d++){
          inputDaily[d] = zmMat[d][a][b];
        }
        FourierCoefs = DFT(inputDaily,nWavesTime);
        outputDaily = IDFT(FourierCoefs);
        for (int d=0; d<dLen; d++){
          outputMat[d][a][b] = outputDaily[d];
        }
      }
    }

    
    
    NcVar *checkDFTThresh = fileout.add_var("THRESHOLD_DFT",ncDouble,outTime,outLat,outLon);
    checkDFTThresh->set_cur(0,0,0);
    checkDFTThresh->put(&(outputMat[0][0][0]),dLen,latLen,lonLen);

    avgin.close();
    fileout.close();
  }
  catch (Exception &e){
    std::cout<<e.ToString()<<std::endl;
  }
}
