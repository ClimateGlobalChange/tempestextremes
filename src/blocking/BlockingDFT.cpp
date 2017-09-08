///////////////////////////////////////////////////
///
///           \file BlockingDFT.cpp
///           \author Marielle Pinheiro
///           \date May 8, 2017
///
//////////////////////////////////////////////////

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
#include <complex>
#include <cstring>

int main(int argc, char ** argv){
  try{
    NcError error(NcError::silent_nonfatal);
    std::string fileList, outFile,avgName,varName,  tname, latname, lonname;
    int nWaves,startday,endday;
    BeginCommandLine()
      CommandLineString(fileList,"inlist","");
      CommandLineString(outFile,"out","");
      CommandLineString(avgName,"avgname","");
      CommandLineString(varName,"varname","");
      CommandLineInt(nWaves, "ncoef",12);
      CommandLineInt(startday,"startday",1);
      CommandLineInt(endday,"endday",365);
      CommandLineString(tname,"tname","time");
      CommandLineString(latname,"latname","lat");
      CommandLineString(lonname,"lonname","lon");
      ParseCommandLine(argc,argv);
    EndCommandLine(argv)
    AnnounceBanner();

    if (fileList == ""){
      _EXCEPTIONT("No file list provided.");
    }
    if (outFile==""){
      _EXCEPTIONT("Need to provide output file name.");
    }
    if (varName==""){
      _EXCEPTIONT("Need to provide input variable name (--varname).");
    }
    if (avgName == ""){
      _EXCEPTIONT("Need to provide output variable name (--avgname).");
    }
    std::vector<std::string> vecFiles;
    GetInputFileList(fileList,vecFiles);
    int nFiles = vecFiles.size();

    //Open up the first file
    NcFile readin(vecFiles[0].c_str());
    if (!readin.is_valid()){
      _EXCEPTION1("Unable to open file %s for reading",vecFiles[0].c_str());
    }
    //axis lengths
    int tLen,latLen,lonLen;
    tLen = readin.get_dim(tname.c_str())->size();
    latLen = readin.get_dim(latname.c_str())->size();
    lonLen = readin.get_dim(lonname.c_str())->size();

    //Initialize storage array: day, lat, lon
    int yearLen = (endday-startday)+1;
    DataMatrix3D<double> storeMat(yearLen,latLen,lonLen);
    //Initialize counts array: day, lat, lon
    DataMatrix3D<double> countsMat(yearLen, latLen, lonLen);
    //Initialize input datan array: lat, lon
    DataMatrix<double> inputData(latLen,lonLen);
    DataVector<double> timeVec(tLen);

    //Close file (then re-open in loop)
    readin.close();

    //Loop through all available files and store values in appropriate day
    for (int x=0; x<nFiles; x++){
      NcFile readin(vecFiles[x].c_str());
      if (!readin.is_valid()){
        _EXCEPTION1("Unable to open file %s for reading",vecFiles[x].c_str());
      }
      std::cout<<"Reading in "<<vecFiles[x].c_str()<<std::endl;
      tLen = readin.get_dim(tname.c_str()) ->size();
      //std::cout<<"Reading in file #"<<x<<", "<<vecFiles[x].c_str()<<std::endl;
      NcVar *inputVar = readin.get_var(varName.c_str());
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
      //Time and calendar
      int dateYear, dateMonth, dateDay, dateHour, dayIndex;
      bool leap;
      //std::cout<<"Storing values."<<std::endl;
      //Store the values and counts in the matrices at the appropriate day index

      //Input data
      for (int t=0; t<tLen; t++){
        inputData.Initialize(latLen, lonLen);
        inputVar->set_cur(t,0,0);
        inputVar->get(&(inputData[0][0]),1,latLen,lonLen);
        ParseTimeDouble(strTimeUnits, strCalendar, timeVec[t],\
          dateYear, dateMonth, dateDay, dateHour);
        leap = checkFileLeap(strTimeUnits, strCalendar, dateYear,\
          dateMonth, dateDay, dateHour, timeVec[t]);
        if (leap == false){
          dayIndex = DayInYear(dateMonth, dateDay)-1;
          //std::cout<<"Day index is "<<dayIndex<<std::endl;
          for (int a=0; a<latLen; a++){
            for (int b=0; b<lonLen; b++){
              storeMat[dayIndex][a][b]+= inputData[a][b];
              countsMat[dayIndex][a][b]+= 1.;
             // if (a==50 && b==20){
             //   std::cout<<"Store: "<<storeMat[dayIndex][a][b]<< " count: "<< countsMat[dayIndex][a][b]<<std::endl;
            //  }
            }
          }
        }
      }
      std::cout<<"Closing file."<<std::endl;
      readin.close();

    }
    //Now average the values by the counts to get the daily average
    for (int d=0; d<yearLen; d++){
      for (int a=0; a<latLen; a++){
        for (int b=0; b<lonLen; b++){
          if (countsMat[d][a][b] < 1){
            storeMat[d][a][b] = 0.;
          }else{
            storeMat[d][a][b] = storeMat[d][a][b]/countsMat[d][a][b];
          }
        }
      }
    }
    //std::cout<<"Divided values to get average."<<std::endl;

    //Output the Fourier transform of the dataset
    std::vector<double> inputDaily(yearLen);
    std::vector<std::complex <double> > FourierCoefs(yearLen);
    std::vector<double> outputDaily(yearLen);
    DataMatrix3D<double> outputMat(yearLen, latLen, lonLen);
    for (int a=0; a<latLen; a++){
      for (int b=0; b<lonLen; b++){
        for (int d=0; d<yearLen; d++){
          inputDaily[d] = storeMat[d][a][b];
          //std::cout<<"daily input for d="<<d<<" is "<<inputDaily[d]<<std::endl;
        }
        FourierCoefs = DFT(inputDaily, nWaves);
        outputDaily = IDFT(FourierCoefs);
        for (int d=0; d<yearLen; d++){
          outputMat[d][a][b] = outputDaily[d];
        }
      }
    }


/*    DataMatrix3D<double> outputMatZonal(yearLen,latLen,lonLen);
    DataMatrix <double> zonalAvgMat(yearLen,latLen);
    for (int d=0; d<yearLen; d++){
      for (int a=0; a<latLen; a++){
        for (int b=0; b<lonLen; b++){
          zonalAvgMat[d][a]+=storeMat[d][a][b]/lonLen;
        }
      }
    }

    //Matrix for the transformed Fourier value
    //Calculate at each lat,lon
    for (int a=0; a<latLen; a++){
      for (int d=0; d<yearLen; d++){
        inputDaily[d] = zonalAvgMat[d][a];
      }
      FourierCoefs = DFT(inputDaily, nWaves);
      outputDaily = IDFT(FourierCoefs);
      for (int d=0; d<yearLen; d++){
        for (int b=0; b<lonLen; b++){
          outputMatZonal[d][a][b] = outputDaily[d];
        }
      }
    }
*/
    //Copy existing values from reference file to output file
    NcFile refFile(vecFiles[0].c_str());
    NcDim *inLatDim = refFile.get_dim(latname.c_str());
    NcDim *inLonDim = refFile.get_dim(lonname.c_str());
    NcVar *inLatVar = refFile.get_var(latname.c_str());
    NcVar *inLonVar = refFile.get_var(lonname.c_str());

    NcFile outfile(outFile.c_str(),NcFile::Replace,NULL,0,NcFile::Offset64Bits);
    NcDim *outTime = outfile.add_dim(tname.c_str(),yearLen);
    DataVector<int> tVals(yearLen);
    for (int t=startday; t<=endday; t++){
      tVals[t] = t;
    }
    NcVar *outTimeVar = outfile.add_var(tname.c_str(),ncInt,outTime);
    outTimeVar->set_cur((long)0);
    outTimeVar->put(&(tVals[0]),yearLen);

    //Add time units
    outTimeVar->add_att("units","days since 0001-01-01");

    NcDim *outLat = outfile.add_dim(latname.c_str(),latLen);
    NcDim *outLon = outfile.add_dim(lonname.c_str(),lonLen);
    NcVar *outLatVar = outfile.add_var(latname.c_str(),ncDouble,outLat);
    NcVar *outLonVar = outfile.add_var(lonname.c_str(),ncDouble,outLon);

    copy_dim_var(inLatVar,outLatVar);
    copy_dim_var(inLonVar,outLonVar);

    NcVar *outAvgVar = outfile.add_var(avgName.c_str(), ncDouble, outTime, outLat, outLon);
    outAvgVar->set_cur(0,0,0);
    outAvgVar->put(&(outputMat[0][0][0]),yearLen,latLen,lonLen);
/*
    std::string zonalAvgName = avgName.append("_NO_SMOOTH");
    NcVar *outZonalVar = outfile.add_var(zonalAvgName.c_str(),ncDouble,outTime,outLat,outLon);
    outZonalVar->set_cur(0,0,0);
    outZonalVar->put(&(storeMat[0][0][0]),yearLen,latLen,lonLen);
*/
    refFile.close();
    outfile.close();

  }
  catch (Exception &e){
    std::cout<<e.ToString()<<std::endl;
  }
}
