////////////////////////////////////
///	\file ExtractTimeStep.cpp
///	\author Marielle Pinheiro
///	\version April 18th, 2019

/*This code extracts the time step of interest for every day
within the file (for example, 12Z).*/

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
#include <string>
#include <cstdio>
#include <vector>
#include <iostream>
#include <sstream>

int main(int argc,char **argv){
	NcError error(NcError::verbose_nonfatal);
	try{
		std::string fileList;
		std::string varName;
		std::string tname,latname,lonname;
		int hourInt;
		
		BeginCommandLine()
			CommandLineString(fileList,"inlist","");
			CommandLineString(varName,"varname","");
			CommandLineString(tname,"tname","time");
			CommandLineString(latname,"latname","lat");
			CommandLineString(lonname,"lonname","lon");
			CommandLineInt(hourInt,"hour",-9999);
			ParseCommandLine(argc,argv);
		EndCommandLine(argv)
		AnnounceBanner();
		if (fileList==""){
			_EXCEPTIONT("Need to provide file list");
		}
		if (std::fabs(hourInt + 9999)<0.01){
			_EXCEPTIONT("Need to provide hour");
		}
		if (varName==""){
			_EXCEPTIONT("Need to provide variable name");
		}
		std::vector<std::string> InputFiles;
		GetInputFileList(fileList,InputFiles);
		int nFiles = InputFiles.size();
		size_t pos,len;
		//Start parsing files
		for (int x=0; x<nFiles; x++){
			//Parse the time variable
			std::cout<<"Reading in file "<<InputFiles[x].c_str()<<std::endl;
			NcFile infile(InputFiles[x].c_str());
			if (!infile.is_valid()){
				_EXCEPTION1("Could not open %s",InputFiles[x].c_str());
			}
			int nTime = infile.get_dim(tname.c_str())->size();
			NcVar * tVar = infile.get_var(tname.c_str());
			if (tVar == NULL){
				_EXCEPTION1("Couldn't find variable %s",tname.c_str());
			}
			DataVector<double> timeVec(nTime);
			tVar->set_cur(long(0));
			tVar->get(&(timeVec[0]),nTime);
			
	//		std::cout<<"Getting time info"<<std::endl;
			//Time units and calendar
			std::string strCalendar;
			std::string strTimeUnits;
			NcAtt * attTime = tVar->get_att("units");
			if (attTime==NULL){
				_EXCEPTIONT("Time variable has no units attribute.");
			}
			strTimeUnits = attTime->as_string(0);
			std::cout<<"Time units are "<<strTimeUnits.c_str()<<std::endl;
			NcAtt * attCal = tVar->get_att("calendar");
	//		std::cout<<"DEBUG: before if statement!"<<std::endl;
			if (attCal==NULL){
				std::cout<<"Calendar attribute does not exist; setting to standard."<<std::endl;
				strCalendar = "standard";
			}else{
				strCalendar = attCal->as_string(0);
			}
	//		std::cout<<"DEBUG: after if statement!"<<std::endl;
			std::cout<<"Calendar is "<<strCalendar.c_str();
			std::vector<int> tIndices;
			int currYear,currMonth,currDay,currHour;
			//Save all time indices with the desired time
			for (int t=0; t<nTime; t++){
				ParseTimeDouble(strTimeUnits,strCalendar,timeVec[t],\
				currYear,currMonth,currDay,currHour);
				if (std::fabs(hourInt-currHour)<0.01){
					tIndices.push_back(t);
				//	std::cout<<"added time index "<<t<<" to vector."<<std::endl;
				}
			}

			//Get lat and lon for other dimensions
			int nLat = infile.get_dim(latname.c_str())->size();
			int nLon = infile.get_dim(lonname.c_str())->size();

			//Create new time variable for output
			int nOutTime = tIndices.size();
			DataVector<double> newTimeVec(nOutTime);
			DataMatrix3D<double> outMat(nOutTime,nLat,nLon);
			//Read in the original variable
			NcVar * inVar = infile.get_var(varName.c_str());
			int currT;
			for (int d=0; d<nOutTime; d++){
				currT = tIndices[d];
				inVar->set_cur(currT,0,0);
				inVar->get(&(outMat[d][0][0]),1,nLat,nLon);
				newTimeVec[d] = timeVec[currT];
			}
			//Write the outfile
			std::string strOutFile;
			std::string suff ="_"+ std::to_string(hourInt) + "Z.nc";
			pos = InputFiles[x].find(".nc");
			len = InputFiles[x].length();
			strOutFile = InputFiles[x].replace(pos,len,suff);
			NcFile outfile(strOutFile.c_str(),NcFile::Replace,NULL,0,NcFile::Offset64Bits);
			std::cout<<"Writing "<<strOutFile.c_str()<<" to file."<<std::endl;
			NcDim * outTime = outfile.add_dim(tname.c_str(),nOutTime);
			NcVar * outTimeVar = outfile.add_var(tname.c_str(),ncDouble,outTime);
			outTimeVar->set_cur(long(0));
			outTimeVar->put(&(newTimeVec[0]),nOutTime);
			outTimeVar->add_att("units",strTimeUnits.c_str());
			outTimeVar->add_att("calendar",strCalendar.c_str());

			NcDim * outLat = outfile.add_dim(latname.c_str(),nLat);
			NcDim * outLon = outfile.add_dim(lonname.c_str(),nLon);
			NcVar * refLat = infile.get_var(latname.c_str());
			NcVar * refLon = infile.get_var(lonname.c_str());
			NcVar * outLatVar = outfile.add_var(latname.c_str(),ncDouble,outLat);
			NcVar * outLonVar = outfile.add_var(lonname.c_str(),ncDouble,outLon);
			copy_dim_var(refLat,outLatVar);
			copy_dim_var(refLon,outLonVar);
			NcVar * outVar = outfile.add_var(varName.c_str(),ncDouble,outTime,outLat,outLon);
			outVar->set_cur(0,0,0);
			outVar->put(&(outMat[0][0][0]),nOutTime,nLat,nLon);
			outfile.close();
			infile.close();
		}

	}
	catch (Exception & e){
		std::cout<<e.ToString()<<std::endl;
	}
}
