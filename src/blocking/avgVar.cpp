//////////////////////////
//
//  \file avgVar.cpp
//
//  \author Marielle Pinheiro
//
//  \date November 15, 2016

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

//Outputs a matrix that is the SUM of the elements along the time axis
void Varto2Dsum(NcVar *inVar, int v, DataMatrix3D<double> & outMat){
    int tLen,latLen,lonLen;
    tLen = inVar->get_dim(0)->size();
    latLen = inVar->get_dim(1)->size();
    lonLen = inVar->get_dim(2)->size();

    DataMatrix3D<double> inMat(tLen,latLen,lonLen);
	//DataMatrix<double> outMat(latLen,lonLen);
    inVar->set_cur(0,0,0);
    inVar->get((&inMat[0][0][0]), tLen,latLen,lonLen);
	
    for (int t=0; t<tLen; t++){
        for (int i=0; i<latLen; i++){
            for (int j=0; j<lonLen; j++){
			    outMat[i][j][v]+=inMat[t][i][j];
            }
        }
    }
}

int main(int argc, char ** argv){
	try{
		std::string fileList, outFile, varList;
		
		BeginCommandLine()
			CommandLineString(fileList,"inlist","");
			CommandLineString(outFile,"out","");
			CommandLineString(varList,"varlist","");
			ParseCommandLine(argc,argv);
		EndCommandLine(argv)
		AnnounceBanner();
		if (fileList == ""){
			_EXCEPTIONT("No file list provided");
		}
		if (outFile == ""){
			_EXCEPTIONT("Need to provide output file name");
		}
		if (varList==""){
			_EXCEPTIONT("Need to provide at least one variable");
		}
		//file list vector
		std::vector<std::string> vecFiles;
		GetInputFileList(fileList,vecFiles);
		
		//variable list vector
		std::vector<std::string> varVec;
	    std::string delim = ",";
	    size_t pos = 0;
	    std::string token;
		
	    while((pos = varList.find(delim)) != std::string::npos){
	      token = varList.substr(0,pos);
	      varVec.push_back(token);
	      varList.erase(0,pos + delim.length());
	    }
	    varVec.push_back(varList);
		
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
		
		//Create your working matrix; 
		// The third axis is the length of your variable vector
		// Use a temporary 2D matrix for function output
		int varLen = varVec.size();
		DataMatrix3D<double> storeMat(latLen,lonLen,varLen);

		//Open the first file and add each variable to working matrix
		for (int v=0; v<varLen; v++){
			NcVar * invar = readin.get_var(varVec[v].c_str());
			Varto2Dsum(invar, v, storeMat);
		}
		
		//Continue adding values to the matrix
		//increment the length of the time axis to accommodate
		for (int f=1; f<vecFiles.size(); f++){
		    NcFile addread(vecFiles[f].c_str());
			for (int v=0; v<varLen; v++){
				NcVar * invar = addread.get_var(varVec[v].c_str());
				tLen += addread.get_dim("time")->size();
				Varto2Dsum(invar, v, storeMat);
			}
			addread.close();	
		}
		double invtLen = 1./((double) tLen);
		//Now divide by the total time length
		for (int a=0; a<latLen; a++){
			for (int b=0; b<lonLen; b++){
				for (int v=0; v<varLen; v++){
					storeMat[a][b][v] = storeMat[a][b][v]*invtLen;
					if (a==5 && b==10){
						std::cout<<" now "<<storeMat[a][b][v]<<std::endl;
					}
				}
			}
		}
		
		//output the variables to file
	    NcFile readout(outFile.c_str(),NcFile::Replace, NULL,0,NcFile::Offset64Bits);
	    NcDim * outLat = readout.add_dim("lat", latLen);
	    NcDim * outLon = readout.add_dim("lon", lonLen);
		NcDim * varvar = readout.add_dim("varnum",varLen);
	    NcVar * outLatVar = readout.add_var("lat",ncDouble,outLat);
	    NcVar * outLonVar = readout.add_var("lon",ncDouble,outLon);
	    std::cout<<"Copying dimension attributes."<<std::endl;
	    copy_dim_var(latVar,outLatVar);
	    copy_dim_var(lonVar,outLonVar);
		
		/*testoutvar = readout.add_var("vars_Test",ncDouble,latLen,lonLen,varLen);
		testoutvar->set_cur(0,0,0);
		testoutvar->put(&(storeMat[0][0][0]),latLen,lonLen,varLen);*/
		//output each of the variables
		DataMatrix<double>putMat(latLen,lonLen);
		for (int v=0; v<varLen; v++){
			NcVar * outvar = readout.add_var(varVec[v].c_str(),ncDouble,outLat,outLon);
			outvar->set_cur(0,0);
			std::cout<<"Value being put in at v="<<v<<" is "<<storeMat[5][10][v]<<std::endl;
			for (int a=0; a<latLen;a++){
				for (int b=0; b<lonLen;b++){
					putMat[a][b]=storeMat[a][b][v];
				}
			}
			outvar->put(&(putMat[0][0]),latLen,lonLen);
		}
		readout.close();
		readin.close();
	}
	catch (Exception &e){
		std::cout<<e.ToString()<<std::endl;
	}
}