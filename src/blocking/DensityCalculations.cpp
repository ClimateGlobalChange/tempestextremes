///////////////////////////////////////////////////
//
//        \file densityCalculations.cpp 
//
//        \author Marielle Pinheiro
//
//        \version 1.0 June 12, 2015

#include "NetCDFUtilities.h"
#include "netcdfcpp.h"
#include "DataVector.h"
#include "DataMatrix3D.h"
#include "DataMatrix.h"
#include "Announce.h"
#include "CommandLine.h"
#include "Exception.h"
#include "BlockingUtilities.h"

#include <cstdlib>
#include <cmath>
#include <cstring>

void densCalc(NcVar * inVar,
              DataMatrix<double> & outMat){

  int tLen,latLen,lonLen;
  double invtLen;
  double selfDiv;
  
  tLen = inVar->get_dim(0)->size();
  invtLen = 1.0/((double) tLen);
  latLen = inVar->get_dim(1)->size();
  lonLen = inVar->get_dim(2)->size();

  DataMatrix3D<double> inMat(tLen,latLen,lonLen);
  inVar->set_cur(0,0,0);
  inVar->get((&inMat[0][0][0]), tLen,latLen,lonLen);

  //DataMatrix<double> outMat(latLen,lonLen);
  for (int i=0; i<latLen; i++){
    for (int j=0; j<lonLen; j++){
      outMat[i][j] = 0.0;
    }
  }

  //Divide all cells by self and average over time
  for (int t=0; t<tLen; t++){
    for (int i=0; i<latLen; i++){
      for (int j=0; j<lonLen; j++){
        if (std::fabs(inMat[t][i][j])>0.0){
          outMat[i][j]+=1.0;
        }
      }
    }
  }

  for (int i=0; i<latLen; i++){
    for (int j=0; j<lonLen; j++){
      outMat[i][j]*=invtLen;
    }
  }
 // outVar->set_cur(0,0);
  //outVar->put((&outMat[0][0]),latLen,lonLen);
}

void yearlyStdDev(DataMatrix3D<double> inMat,
                 int nTime,
                 int nLat,
                 int nLon,
                 DataMatrix<double> & outMat){
  DataMatrix<double> meanval(nLat,nLon);
  double invt = 1./((double)nTime);

  for (int t=0; t<nTime; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        meanval[a][b]+=(inMat[t][a][b]*invt);
      }
    }
  }
  double dev;
  double variance;
  for (int a=0; a<nLat; a++){
    for (int b=0; b<nLon; b++){
      dev=0;
      variance=0;
      for (int t=0; t<nTime; t++){
        dev=inMat[t][a][b]-meanval[a][b];
        variance+=(dev*dev*invt);
      }
      outMat[a][b]=std::sqrt(variance);
    }
  }
}

int main(int argc, char ** argv){
  NcError error(NcError::verbose_nonfatal);

  try{
    std::string inFile;
    std::string outFile;
    std::string varName;
    std::string inList;
//    bool calcStdDev;

    BeginCommandLine()
      CommandLineString(inFile, "in", "");
      CommandLineString(inList, "inlist","");
      CommandLineString(varName, "var", "");
      CommandLineString(outFile, "out", "");
  //    CommandLineBool(calcStdDev, "std");
      ParseCommandLine(argc, argv);
    EndCommandLine(argv)
    AnnounceBanner();
  

    if ((inFile != "") && (inList != "")){
      _EXCEPTIONT("Can only open one file (--in) or list (--inlist).");
    }

    //file list vector
    std::vector<std::string> vecFiles;
    if (inFile != ""){
      vecFiles.push_back(inFile);
    }
    if (inList != ""){
      GetInputFileList(inList,vecFiles);
    }

    //open up first file
    NcFile readin(vecFiles[0].c_str());
    if (!readin.is_valid()){
      _EXCEPTION1("Unable to open file %s for reading",\
        vecFiles[0].c_str());
    }
    int tLen,latLen,lonLen;

    NcDim * time = readin.get_dim("time");
    tLen = time->size();
    NcVar * timeVar = readin.get_var("time");

    NcDim * lat = readin.get_dim("lat");
    latLen = lat->size();
    NcVar * latVar = readin.get_var("lat");

    NcDim * lon = readin.get_dim("lon");
    lonLen = lon->size();
    NcVar * lonVar = readin.get_var("lon");

    //read input variable
    NcVar * inVar = readin.get_var(varName.c_str());
    //Create output matrix
    DataMatrix<double> outMat(latLen,lonLen);
    densCalc(inVar,outMat);
    //Option for calculating the yearly standard deviation 
 /*   if (calcStdDev){
      for (int a=0; a<latLen; a++){
        for (int b=0; b<lonLen; b++){
          storeMat[0][a][b] = outMat[a][b];
        }
      }
    }
*/
    //If multiple files, add these values to the output
    if (vecFiles.size()>1){
      DataMatrix<double> addMat(latLen,lonLen);
      std::cout<<"There are "<<vecFiles.size()<<" files."<<std::endl;
      for (int v=1; v<vecFiles.size(); v++){
        NcFile addread(vecFiles[v].c_str());
        NcVar * inVar = addread.get_var(varName.c_str());
        densCalc(inVar,addMat);
        for (int a=0; a<latLen; a++){
          for (int b=0; b<lonLen; b++){
            outMat[a][b]+=addMat[a][b];
          }
        }
/*        if (calcStdDev){
          for (int a=0; a<latLen; a++){
            for (int b=0; b<lonLen; b++){
              storeMat[v][a][b] = addMat[a][b];
            }
          }
        }*/
        addread.close(); 
      }
      //Divide output by number of files
      double div = 1./((double) vecFiles.size());
      for (int a=0; a<latLen; a++){
        for (int b=0; b<lonLen; b++){
          outMat[a][b]*=div;
        }
      }
    }
 
    NcFile readout(outFile.c_str(),NcFile::Replace, NULL,0,NcFile::Offset64Bits);
    NcDim * outLat = readout.add_dim("lat", latLen);
    NcDim * outLon = readout.add_dim("lon", lonLen);
    NcVar * outLatVar = readout.add_var("lat",ncDouble,outLat);
    NcVar * outLonVar = readout.add_var("lon",ncDouble,outLon);
    std::cout<<"Copying dimension attributes."<<std::endl;
    copy_dim_var(latVar,outLatVar);
    copy_dim_var(lonVar,outLonVar);

    std::cout<<"Creating density variable."<<std::endl;

    NcVar * densVar = readout.add_var("dens",ncDouble,outLat,outLon);
    densVar->set_cur(0,0);
    densVar->put((&outMat[0][0]),latLen,lonLen);

/*    if (calcStdDev){
      NcVar * stdDevVar = readout.add_var("stddev", ncDouble,outLat,outLon);
      DataMatrix<double> stdDevMat(latLen,lonLen);
      yearlyStdDev(storeMat,vecFiles.size(),latLen,lonLen,stdDevMat);
      stdDevVar->set_cur(0,0);
      stdDevVar->put(&(stdDevMat[0][0]),latLen,lonLen);
      std::cout<<" created sd variable"<<std::endl;

    }
*/
    readout.close();
    readin.close();
  } 
  catch (Exception &e){
    std::cout<<e.ToString()<<std::endl;
  }
}
