////////////////////////////////////
///     \file DailyAvg.cpp
///     \author Marielle Pinheiro
///     \version April 22nd, 2019

/*This code averages values for each day in the file*/

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

int main(int argc,char**argv){
	NcError error(NcError::verbose_nonfatal);
	try{
                std::string fileList;
                std::string varName;
                std::string tname,latname,lonname,levname;
		bool is4d,isHpa,g2z;
                BeginCommandLine()
                        CommandLineString(fileList,"inlist","");
                        CommandLineString(varName,"varname","");
                        CommandLineString(tname,"tname","time");
			CommandLineString(levname,"levname","lev");
                        CommandLineString(latname,"latname","lat");
                        CommandLineString(lonname,"lonname","lon");
			CommandLineBool(is4d,"is4d");
			CommandLineBool(isHpa,"hpa");
			CommandLineBool(g2z,"g2z");
                        ParseCommandLine(argc,argv);
                EndCommandLine(argv)
                AnnounceBanner();
                if (fileList==""){
                        _EXCEPTIONT("Need to provide file list");
                }
                if (varName==""){
                        _EXCEPTIONT("Need to provide variable name");
                }
                std::vector<std::string> InputFiles;
                GetInputFileList(fileList,InputFiles);
                int nFiles = InputFiles.size();
                size_t pos,len;

		double gmult = 1.;
		if (g2z){
			gmult = 1./9.8;
		}
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

                       // std::cout<<"Getting time info"<<std::endl;
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
                        //std::cout<<"DEBUG: before if statement!"<<std::endl;
                        if (attCal==NULL){
                                std::cout<<"Calendar attribute does not exist; setting to standard."<<std::endl;
                                strCalendar = "standard";
                        }else{
                                strCalendar = attCal->as_string(0);
                        }
                        //std::cout<<"DEBUG: after if statement!"<<std::endl;
                        std::cout<<"Calendar is "<<strCalendar.c_str();
                        std::vector<int> tIndices;
                        int prevYear=-9999;
			int prevMonth=-9999;
			int prevDay=-9999;
			int prevHour=-9999;
			int currYear,currMonth,currDay,currHour;
			//Vector of time step with each new day
			for (int t=0; t<nTime; t++){
				ParseTimeDouble(strTimeUnits,strCalendar,timeVec[t],\
				currYear,currMonth,currDay,currHour);
		//		std::cout<<"prevdate is "<<prevYear<<prevMonth<<prevDay<<" and currdate is "<<currYear<<currMonth<<currDay<<std::endl;
				if (std::fabs(currDay-prevDay)>0.01){
					tIndices.push_back(t);
					std::cout<<"Added time step "<<t<<" ("<<currYear<<currMonth<<currDay<<") to vector"<<std::endl;
					prevYear = currYear;
					prevMonth=currMonth;
					prevDay = currDay;
					prevHour = currHour;
				}
			}
                        //Get lat and lon for other dimensions
                        int nLat = infile.get_dim(latname.c_str())->size();
                        int nLon = infile.get_dim(lonname.c_str())->size();

                        //Create new time variable for output
                        int nOutTime = tIndices.size();
			std::cout<<"Output will be "<<nOutTime<<" in time dimension"<<std::endl;
                        DataVector<double> newTimeVec(nOutTime);
			int pIndex = 10000000;
    			if (is4d){
      				int nLev = infile.get_dim(levname.c_str())->size();
    				NcVar *levvar = infile.get_var(levname.c_str());
      				DataVector<double> pVec(nLev);
      				levvar->set_cur(long(0));
      				levvar->get(&pVec[0],nLev);

      				double pval = 50000.0;
      				if (isHpa){
        				pval=500.;
      				}
      				for (int p=0; p<nLev; p++){
        				if (std::fabs(pVec[p]-pval)<0.0001){
          					pIndex = p;
                				break;
              				}
      				}
      				if (pIndex > 999999){
        				_EXCEPTIONT("Could not identify correct pressure level. Check file.");
    				}
		    	}

                        //Read in the original variable
                        NcVar * inVar = infile.get_var(varName.c_str());
			DataMatrix3D<double> datStore(nOutTime,nLat,nLon);
			DataMatrix3D<double> countStore(nOutTime,nLat,nLon);
			DataMatrix<double> currSlice(nLat,nLon);
			//Add each time slice to the appropriate day
			int ds,de;
			for (int d=0; d<nOutTime; d++){
				ds=tIndices[d];
				newTimeVec[d]= timeVec[ds];
				if (d<(nOutTime-1)){
					de=tIndices[d+1]-1;
				}else{
					de=nTime-1;
				}
				std::cout<<"ds is "<<ds<<" and de is "<<de<<std::endl;
				for (int s=ds; s<=de; s++){
				//	std::cout<<"s is currently "<<s<<std::endl;
					//Set the time slice value
					if (is4d){
						inVar->set_cur(s,pIndex,0,0);
						inVar->get(&(currSlice[0][0]),1,1,nLat,nLon);
					}else{
						inVar->set_cur(s,0,0);
						inVar->get(&(currSlice[0][0]),1,nLat,nLon);
					}
					for (int a=0; a<nLat; a++){
						for (int b=0; b<nLon; b++){
							datStore[d][a][b]+=currSlice[a][b]*gmult;
							countStore[d][a][b]+=1.;
						}
					}
				}
			}
			DataMatrix3D<double> avgStore(nOutTime,nLat,nLon);
			//Now average all the values
			for (int d=0; d<nOutTime; d++){
				for (int a=0; a<nLat; a++){
					for (int b=0; b<nLon; b++){
						avgStore[d][a][b] = datStore[d][a][b]/countStore[d][a][b];
					}
				}
			}

			//Write the outfile
                        std::string strOutFile;
                        std::string suff ="_dailyavg.nc";
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
                        outVar->put(&(avgStore[0][0][0]),nOutTime,nLat,nLon);
                        outfile.close();
                        infile.close();


		}
	}
	catch (Exception & e){
		std::cout<<e.ToString()<<std::endl;
	}
}
