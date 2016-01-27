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
#include "blockingUtilities.h"

#include <cstdlib>
#include <cmath>
#include <cstring>

void densCalc(NcVar * inVar,
              NcVar * outVar){

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

  DataMatrix<double> outMat(latLen,lonLen);
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
  outVar->set_cur(0,0);
  outVar->put((&outMat[0][0]),latLen,lonLen);
}

int main(int argc, char ** argv){
  NcError error(NcError::verbose_nonfatal);

  try{
    std::string inFile;
    std::string outFile;
    std::string varName;

    BeginCommandLine()
      CommandLineString(inFile, "in", "");
      CommandLineString(varName, "var", "");
      CommandLineString(outFile, "out", "");
      ParseCommandLine(argc, argv);
    EndCommandLine(argv)
    AnnounceBanner();
  
    NcFile readin(inFile.c_str());
    if (!readin.is_valid()){
      _EXCEPTION1("Unable to open file %s for reading",\
     inFile.c_str());
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

    NcVar * inVar = readin.get_var(varName.c_str());
    NcFile readout(outFile.c_str(),NcFile::Replace, NULL,0,NcFile::Offset64Bits);
    NcDim * outLat = readout.add_dim("lat", latLen);
    NcDim * outLon = readout.add_dim("lon", lonLen);
    NcVar * outLatVar = readout.add_var("lat",ncDouble,outLat);
    NcVar * outLonVar = readout.add_var("lon",ncDouble,outLon);
    copy_dim_var(latVar,outLatVar);
    copy_dim_var(lonVar,outLonVar);

    NcVar * densVar = readout.add_var("dens",ncDouble,outLat,outLon);
    densCalc(inVar,densVar);

    readout.close();
    readin.close();
  } 
  catch (Exception &e){
    std::cout<<e.ToString()<<std::endl;
  }
}
