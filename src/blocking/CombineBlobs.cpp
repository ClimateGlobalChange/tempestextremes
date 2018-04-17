/////////////////////////////////////
///     \file CombineBlobs.cpp
///     \author Kyle Stachowicz
///     \version May 23, 2017

#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include "BlockingUtilities.h"
#include "CommandLine.h"
#include "DataMatrix3D.h"
#include "netcdfcpp.h"

int main(int argc, char** argv) {
  std::string inListFileName;
  std::string outFileName;

  std::string timeDimName, latDimName, lonDimName;

  std::string blobVarName;

  BeginCommandLine() {
    CommandLineString(inListFileName, "inlist", "");
    CommandLineString(outFileName, "out", "");

    CommandLineString(timeDimName, "timeDim", "time");
    CommandLineString(latDimName, "latDim", "lat");
    CommandLineString(lonDimName, "lonDim", "lon");

    CommandLineString(blobVarName, "blobVar", "");
    ParseCommandLine(argc, argv);
  }
  EndCommandLine(argv);
  AnnounceBanner();

  if (blobVarName.length() == 0) {
    std::cerr << "Error: no blob variable name provided!" << std::endl;
    std::exit(-1);
  }

  if (inListFileName.length() == 0) {
    std::cerr << "Error: no input list provided!" << std::endl;
    std::exit(-1);
  }

  std::vector<std::string> inFileNames;
  GetInputFileList(inListFileName,inFileNames);
/*  {
    std::ifstream inListFile{inListFileName};
    std::string fileNameBuf;
    while (std::getline(inListFile, fileNameBuf)) {
      inFileNames.push_back(fileNameBuf);
    }
  }
*/

  int nFiles = inFileNames.size();

/*
  std::vector<std::unique_ptr<NcFile>> inFiles;
  std::string fileName;
//  for (auto& fileName : inFileNames) {
  for (int x=0; x<nFiles; x++){
    fileName = inFileNames[x];
    try {
      inFiles.emplace_back(new NcFile(fileName.c_str()));
    } catch (...) {
      std::cerr << fileName << " failed to load!" << std::endl;
    }
  }

  if (outFileName.length() == 0) {
    std::cerr << "Error: no output filename provided!" << std::endl;
    std::exit(-1);
  }
*/

  NcFile refFile(inFileNames[0].c_str());
  size_t timeDimSize = refFile.get_dim(timeDimName.c_str())->size();
  size_t latDimSize = refFile.get_dim(latDimName.c_str())->size();
  size_t lonDimSize = refFile.get_dim(lonDimName.c_str())->size();

  NcFile outFile{outFileName.c_str(), NcFile::Replace};
  NcDim* outTimeDim = outFile.add_dim(timeDimName.c_str(), timeDimSize);
  NcDim* outLatDim = outFile.add_dim(latDimName.c_str(), latDimSize);
  NcDim* outLonDim = outFile.add_dim(lonDimName.c_str(), lonDimSize);

  {
    NcVar* outTimeVar =
        outFile.add_var(timeDimName.c_str(), ncDouble, outTimeDim);
    NcVar* outLatVar = outFile.add_var(latDimName.c_str(), ncDouble, outLatDim);
    NcVar* outLonVar = outFile.add_var(lonDimName.c_str(), ncDouble, outLonDim);

    NcVar* inTimeVar = refFile.get_var(timeDimName.c_str());
    NcVar* inLatVar = refFile.get_var(latDimName.c_str());
    NcVar* inLonVar = refFile.get_var(lonDimName.c_str());

    copy_dim_var(inTimeVar, outTimeVar);
    copy_dim_var(inLatVar, outLatVar);
    copy_dim_var(inLonVar, outLonVar);
  }

  refFile.close();
  NcVar* outBlobVar = outFile.add_var(blobVarName.c_str(), ncInt, outTimeDim,
                                      outLatDim, outLonDim);

  DataMatrix3D<int> outData(timeDimSize, latDimSize, lonDimSize);
  DataMatrix3D<int> inData(timeDimSize, latDimSize, lonDimSize);

//  for (auto& inFile : inFiles) {
  for (int x=0; x<nFiles; x++){
    NcFile inFile(inFileNames[x].c_str());
    NcDim* timeDim = inFile.get_dim(timeDimName.c_str());
    assert(timeDimSize == timeDim->size());

    NcDim* latDim = inFile.get_dim(latDimName.c_str());
    assert(latDimSize == latDim->size());

    NcDim* lonDim = inFile.get_dim(lonDimName.c_str());
    assert(lonDimSize == lonDim->size());

    NcVar* blobVar = inFile.get_var(blobVarName.c_str());
    blobVar->set_cur(0, 0, 0);
    blobVar->get(&(inData[0][0][0]), timeDimSize, latDimSize, lonDimSize);

    for (size_t t = 0; t < timeDimSize; t++) {
      for (size_t lat = 0; lat < latDimSize; lat++) {
        for (size_t lon = 0; lon < lonDimSize; lon++) {
          outData[t][lat][lon] = outData[t][lat][lon] || inData[t][lat][lon];
        }
      }
    }
    inFile.close();
  }

  outBlobVar->put(&(outData[0][0][0]), timeDimSize, latDimSize, lonDimSize);

  outFile.close();
}
