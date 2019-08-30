/////////////////////////////////////
///     \file DetrendHeights.cpp
///     \author Marielle Pinheiro
///     \version March 25, 2019


/*This code takes the linear regression 
coefficients from the calcLinReg.py output
and detrends Z500 data for each day in the
year, with the values at the specified
central year acting as the new 0 axis for 
the detrended data.
*/

#include "BlockingUtilities.h"
#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "TimeObj.h"

#include "DataArray1D.h"
#include "DataArray2D.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <sstream>

int main(int argc, char **argv){
    NcError error(NcError::verbose_nonfatal);
    try{
        std::string fileList,linfile;
        std::string varName;
        std::string tname, latname, lonname, levname;
        bool is4d, isHpa, GHtoZ;
        int centerYear;
        //,endYear,centerYear;
        BeginCommandLine()
            CommandLineString(fileList,"inlist","");
            CommandLineString(linfile,"linfile","");
            CommandLineString(varName,"varname","");
            //CommandLineInt(startYear,"startyear",-9999);
            //CommandLineInt(endYear,"endyear",-9999);
            CommandLineInt(centerYear,"centeryear",-9999);
            CommandLineString(tname,"tname","time");
            CommandLineString(latname,"latname","lat");
            CommandLineString(lonname,"lonname","lon");
            CommandLineString(levname,"levname","lev");
            CommandLineBool(is4d,"is4D");
            CommandLineBool(isHpa,"hpa");
            CommandLineBool(GHtoZ,"g2z");
        ParseCommandLine(argc, argv);
        EndCommandLine(argv)
        AnnounceBanner();
        if (fileList == ""){
            _EXCEPTIONT("No file list (--inlist) provided");
        }
        if (varName == ""){
            _EXCEPTIONT("No variable name (--varname) specified");
        }

        if (std::fabs(centerYear+9999)<0.001){
            _EXCEPTIONT("Need to provide center year of trend line (--centeryear)");
        }

        //Generate the file list
        std::vector<std::string> InputFiles;
        GetInputFileList(fileList, InputFiles);
        int nFiles = InputFiles.size();   
        //int nYearsInFiles = endYear-startYear+1;
        /*if (std::fabs(centerYear+9999)<0.001){
            if (std::fabs(startYear-9999)<0.001){
                _EXCEPTIONT("Need to provide first year of input files (--startyear)");
            }
            if (std::fabs(endYear-9999)<0.001){
                _EXCEPTIONT("Need to provide last year of input files (--endyear)");
            }
	    std::cout<<"Start year is "<<startYear<<" and end year is "<<endYear<<std::endl;
            //How many years are in the file list?
            int nYearsInFiles = endYear-startYear +1;
            //What is the value of half the years +1?
            centerYear = std::floor(nYearsInFiles/2)+startYear;
            std::cout<<"For "<<nYearsInFiles<<" years, the center year is "<<centerYear<<std::endl;
        }*/
        //Multiplier for GH info
        double ghMult=1.;
        if (GHtoZ){
            ghMult = 1./9.8;
        }
        int pIndex = 1000000;
        if (is4d){
            //Find the pressure level for Z500
            NcFile reffile(InputFiles[0].c_str());
            int nLev = reffile.get_dim(levname.c_str())->size();
            NcVar * levvar = reffile.get_var(levname.c_str());
            DataArray1D<double> pVec(nLev);
            levvar->set_cur(long(0));
            levvar->get(&(pVec[0]),nLev);

            double pval = 50000.;
            if (isHpa){
                pval = 500.;
            }
            for (int x=0; x<nLev; x++){
                if (std::fabs(pVec[x]-pval)<0.001){
                    pIndex=x;
                    break;
                }
            }
            if (pIndex > 99999){
                _EXCEPTIONT("Could not identify correct pressure level. Check file.");
            }

        }
        //Read in the linear regression data
        NcFile linreg(linfile.c_str());
        std::cout<<"Reading in the linear regression data from "<<linfile.c_str();
        int lintime = linreg.get_dim(tname.c_str())->size();
        int nLat = linreg.get_dim(latname.c_str())->size();
        int nLon = linreg.get_dim(lonname.c_str())->size();
        //Time variable for linear regression
        //Note: in units of "days since 0001-01-01"
        //We already know the day in year, which will be used to match the 
        //regression value to the date
        NcVar * lintimevar = linreg.get_var(tname.c_str());
        DataArray1D<double> linTimeVec(lintime);
        lintimevar->set_cur((long)0);
        lintimevar->get(&(linTimeVec[0]),lintime);

        std::string linUnits;
        std::string linCalendar;
        NcAtt * timeUnits = lintimevar->get_att("units");
        if (timeUnits==NULL){
            _EXCEPTIONT("Could not find time units.");
        }
        linUnits = timeUnits->as_string(0);

        NcAtt * timeCal = lintimevar->get_att("calendar");
        linCalendar = timeCal->as_string(0);
        if (timeCal==NULL || strncmp(linCalendar.c_str(),"gregorian",9)==0){
            linCalendar="standard";
        }

        //Open the linear trend variables
        NcVar * slopeVar = linreg.get_var("slope");
        //NcVar * interceptVar = linreg.get_var("intercept");

   

        int nTime;
        std::string strTimeUnits,strCalendar;
        int dateYear,dateMonth,dateDay,dateHour;
        int dayIndex;
        int yearDiff;
        double detrendVal,detrendMean;
        DataArray2D<double> slopeStore(nLat,nLon);
        DataArray2D<double>interceptStore(nLat,nLon);
        DataArray2D<double>varSlice(nLat,nLon);
        size_t pos,len;
        //Begin reading in the files
        for (int x=0; x<nFiles; x++){
            NcFile infile(InputFiles[x].c_str());
            if (!infile.is_valid()){
                _EXCEPTION1("Could not open file %s",InputFiles[x].c_str());
            }
            std::cout<<"Opened "<<InputFiles[x].c_str()<<std::endl;
            NcVar * latvar = infile.get_var(latname.c_str());
            NcVar * lonvar = infile.get_var(lonname.c_str());
            //Time variable
            NcVar * timevar = infile.get_var(tname.c_str());
            nTime = infile.get_dim(tname.c_str())->size();
            DataArray1D<double>timeVec(nTime);
            timevar->set_cur(long(0));
            timevar->get(&(timeVec[0]),nTime);
            //Time units
            NcAtt *attTime = timevar->get_att("units");
            if (attTime == NULL){
                _EXCEPTIONT("Time variable has no units attribute.");
            }
            strTimeUnits = attTime->as_string(0);
            NcAtt *attCal = timevar->get_att("calendar");
	    if (attCal==NULL){
		    strCalendar = "standard";
		 }else{
            strCalendar = attCal->as_string(0);
		}
	    //Read in the height data
            NcVar * heightData = infile.get_var(varName.c_str());
            if (heightData==NULL){
                _EXCEPTIONT("Couldn't read in variable.");
            }
            DataArray3D<double> detrendStore(nTime,nLat,nLon);
            for (int t=0; t<nTime; t++){
                //Read in the time info
                ParseTimeDouble(strTimeUnits,strCalendar,timeVec[t],\
                dateYear,dateMonth,dateDay,dateHour);
                //Get the day of the year for accessing the slope/intercept
                dayIndex = DayInYear(dateMonth,dateDay,strCalendar)-1;
                slopeVar->set_cur(dayIndex,0,0);
                slopeVar->get(&(slopeStore[0][0]),1,nLat,nLon);
                //Difference between current year and starting year of trendline
                yearDiff=dateYear-centerYear;
		std::cout<<"year diff is "<<yearDiff<<std::endl;
                //Get the time slice for the original variable
                if (is4d){
                    heightData->set_cur(t,pIndex,0,0);
                    heightData->get(&(varSlice[0][0]),1,1,nLat,nLon);
                }else{
                    heightData->set_cur(t,0,0);
                    heightData->get(&(varSlice[0][0]),1,nLat,nLon);
                }
                for (int a=0; a<nLat; a++){
                    for (int b=0; b<nLon; b++){
                        //Detrended value is:
                        //D=Z - L 
                        //Where t=0 at the center year and the year difference
                        //is calculated with respect to the center year
                        detrendVal = slopeStore[a][b]* double(yearDiff);
			if (a==10 && b==50){
				std::cout<<"slope is "<<slopeStore[a][b]<<" and detrend val is "<<detrendVal<<std::endl;
			}
                        detrendStore[t][a][b]= varSlice[a][b]*ghMult - detrendVal;
                    }
                }
            }
            //Now write the detrended data to file
            std::string strOutFile;
            pos=InputFiles[x].find(".nc");
            len=InputFiles[x].length();
            strOutFile = InputFiles[x].replace(pos,len,"_detrend_z500.nc");
            NcFile outfile(strOutFile.c_str(),NcFile::Replace,NULL,0,NcFile::Offset64Bits);
            NcDim * outTime = outfile.add_dim(tname.c_str(),nTime);
            NcDim * outLat = outfile.add_dim(latname.c_str(),nLat);
            NcDim * outLon = outfile.add_dim(lonname.c_str(),nLon);

            NcVar * outTimeVar = outfile.add_var(tname.c_str(),ncDouble,outTime);
            NcVar * outLatVar = outfile.add_var(latname.c_str(),ncDouble,outLat);
            NcVar * outLonVar = outfile.add_var(lonname.c_str(),ncDouble,outLon);
            copy_dim_var(timevar,outTimeVar);
            copy_dim_var(latvar,outLatVar);
            copy_dim_var(lonvar,outLonVar);

            NcVar * outVar = outfile.add_var(varName.c_str(),ncDouble,outTime,outLat,outLon);
            outVar->set_cur(0,0,0);
            outVar->put(&(detrendStore[0][0][0]),nTime,nLat,nLon);
            outfile.close();

        }
    }
    catch (Exception & e){
        std::cout<<e.ToString()<<std::endl;
    }
}
