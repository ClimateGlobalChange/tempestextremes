/////////////////////////////////////////////////////////
///
///           \file devIPV.cpp
///           \author Marielle Pinheiro
///           \version April 9, 2015



#include "CLIVAR_block_utilities.h"
#include "CommandLine.h"
#include "Announce.h"
#include "Exception.h"
#include "NetCDFUtilities.h"
#include "netcdfcpp.h"
#include "DataVector.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>

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
void calcDevs(NcVar *inIPV, 
              NcVar *outDev, 
              NcVar *outADev,
              NcVar *outPosIntDev,
              NcVar *avgIPV, 
              NcVar *inTime, 
              NcVar *avgTime,
              NcVar *lat,
              double PVAnom){

  int nTime = inIPV->get_dim(0)->size();
  int nLat = inIPV->get_dim(1)->size();
  int nLon = inIPV->get_dim(2)->size();

//keep sign consistent!
  if (PVAnom<0){
    PVAnom = -PVAnom;
  }

//input IPV
  DataMatrix3D<double> IPVMat(nTime,nLat,nLon);
  inIPV->set_cur(0,0,0);
  inIPV->get(&(IPVMat[0][0][0]),nTime,nLat,nLon);

  DataVector<double> timeVec(nTime);
  inTime->set_cur((long) 0);
  inTime->get(&(timeVec[0]),nTime);

  double tRes = timeVec[1]-timeVec[0];
  int nSteps = 1/tRes;

  std::cout<<"Time resolution is "<<tRes<< "and steps per day is "<<nSteps<<std::endl;
 
//avg IPV
  int avgDay = avgIPV->get_dim(0)->size();
  DataMatrix3D<double> avgMat(avgDay,nLat,nLon);
  avgIPV->set_cur(0,0,0);
  avgIPV->get(&(avgMat[0][0][0]),avgDay,nLat,nLon);

  DataVector<int> avgDayVec(avgDay);
  avgTime->set_cur((long) 0);
  avgTime->get(&(avgDayVec[0]),avgDay);


//Matrix for output data
  DataMatrix3D<double> devMat(nTime,nLat,nLon);
  DataMatrix3D<double> aDevMat(nTime,nLat,nLon);
//Each day needs to correspond to one of the averaged days
//Find the start date that corresponds to one of those days
  double startTime = std::fmod(timeVec[0],365);
  std::cout<<"First day in time axis is day "<<startTime<<std::endl;

//First index in average array
  int avgIndex = int(startTime)-1;

//Number of days in IPV
  int nDays = nTime*tRes;

  std::cout<<"There are "<<nDays<<" days in file."<<std::endl;


//Calculate deviations
  for (int d=0; d<nDays; d++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        for (int t=0; t<nSteps; t++){
          devMat[d*nSteps+t][a][b] = IPVMat[d*nSteps+t][a][b]-avgMat[avgIndex][a][b];
        //  if (a<5 && b<5){
        //    std::cout<< "Dev matrix index is "<<d*nSteps+t<<" and avg index is "\
        //      <<avgIndex<<std::endl;
        //  }
        }
      }
    }
    avgIndex+=1;
  }

  //implement 2-day
  for (int t=0; t<nTime; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        if (t<2*nSteps){
          aDevMat[t][a][b] = devMat[t][a][b];
        }
        else{
          for (int n=0; n<2*nSteps; n++){
            aDevMat[t][a][b]+=devMat[t-n][a][b];
          }
          aDevMat[t][a][b] = aDevMat[t][a][b]/float(2*nSteps);
        }
      }
    }
  }
  outDev->set_cur(0,0,0);
  outDev->put(&(devMat[0][0][0]),nTime,nLat,nLon);

  outADev->set_cur(0,0,0);
  outADev->put(&(aDevMat[0][0][0]),nTime,nLat,nLon);

//Divide matrix by PV anomaly value 
//We are looking for negative anomalies in NH and positive anomalies in SH
  DataVector<double> latVec(nLat);
  lat->set_cur((long) 0);
  lat->get(&(latVec[0]),nLat);

  DataMatrix3D<int> posIntDevs(nTime,nLat,nLon);

  for (int t=0; t<nTime; t++){
    for (int a=0; a<nLat; a++){
      for (int b=0; b<nLon; b++){
        double divDev = aDevMat[t][a][b]/PVAnom;
        //SH: positive anomalies
        if (latVec[a]<0){
          double pos = (divDev+std::fabs(divDev))/2.0;
          posIntDevs[t][a][b] = int(pos);
        }
        //NH: negative anomalies
        else if (latVec[a]>=0){
          double neg = (divDev-std::fabs(divDev))/2.0;
          posIntDevs[t][a][b] = int(neg);
        }
      }
    }
  }
  outPosIntDev->set_cur(0,0,0);
  outPosIntDev->put(&(posIntDevs[0][0][0]),nTime,nLat,nLon);
}





int main(int argc, char **argv){
  NcError error(NcError::verbose_nonfatal);

  try{

    //list of input files for which to calculate deviations
    std::string fileList;
    //file that holds averages
    std::string avgName;

    BeginCommandLine()
      CommandLineString(fileList, "in", "");
      CommandLineString(avgName, "avg", "");
      ParseCommandLine(argc, argv);
    EndCommandLine(argv);
    AnnounceBanner();

   //Create list of input files
    std::vector<std::string> InputFiles;
    GetInputFileList(fileList, InputFiles);
    int nFiles = InputFiles.size();

    //Open averages file
    NcFile avgFile(avgName.c_str());

    //time data (average)
    int avgTime = avgFile.get_dim("time")->size();
    NcVar *avgTimeVals = avgFile.get_var("time");
    
    //averaged IPV
    NcVar *AIPVData = avgFile.get_var("AIPV");
    double anomVal = std::pow(1.3,-6);
    //Open IPV files

    for (int x=0; x<nFiles; x++){

      NcFile infile(InputFiles[0].c_str());
      NcDim *tDim = infile.get_dim("time");
      NcDim *latDim = infile.get_dim("lat");
      NcDim *lonDim = infile.get_dim("lon");
      
      int nTime = tDim->size();
      int nLat = latDim->size();
      int nLon = lonDim->size();

      NcVar *inTime = infile.get_var("time");
      NcVar *inLat = infile.get_var("lat");
      NcVar *inLon = infile.get_var("lon");
      NcVar *IPVData = infile.get_var("IPV");

      //Create output file that corresponds to IPV data
      std::string strOutFile = InputFiles[x].replace(InputFiles[x].end()-3,\
        InputFiles[x].end(), "_devs.nc");

      NcFile outfile(strOutFile.c_str(), NcFile::Replace, NULL,0,NcFile::Offset64Bits);
      NcDim *tDimOut = outfile.add_dim("time",nTime);
      NcDim *latDimOut = outfile.add_dim("lat",nLat);
      NcDim *lonDimOut = outfile.add_dim("lon",nLon);

      NcVar *tVarOut = outfile.add_var("time",ncDouble,tDimOut);
      copy_dim_var(inTime,tVarOut);
      NcVar *latVarOut = outfile.add_var("lat",ncDouble,latDimOut);
      copy_dim_var(inLat,latVarOut);
      NcVar *lonVarOut = outfile.add_var("lon",ncDouble,lonDimOut);
      copy_dim_var(inLon,lonVarOut);

      //Create variables for Deviations
      NcVar *devOut = outfile.add_var("DIPV",ncDouble,tDimOut,latDimOut,lonDimOut);
      NcVar *aDevOut = outfile.add_var("ADIPV",ncDouble,tDimOut,latDimOut,lonDimOut);
      NcVar *devIntOut = outfile.add_var("INT_DIPV",ncInt,tDimOut,latDimOut,lonDimOut);

      calcDevs(IPVData, devOut, aDevOut, devIntOut, AIPVData, inTime, avgTimeVals, inLat, anomVal);
 
      infile.close();
      outfile.close();
    }

  }catch (Exception & e){
    std::cout<<e.ToString() <<std::endl;
  }

}
