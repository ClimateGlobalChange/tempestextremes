/////////////////////////////////////////////////////////
///
///           \file blockingDevs.cpp
///           \author Marielle Pinheiro
///           \version June 1, 2015

/*This code is the third step in the potential vorticity code
based on the Schwierz et al 2004 paper. It calculates deviations 
from the average produced in the previous step (DIPV) and filters them with
2-day smoothing (ADIPV), then divides ADIPV by the specified deviation 
and outputs integer values (INT_ADIPV) which can then be used by StitchBlobs

*/

#include "BlockingUtilities.h"
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


int main(int argc, char **argv){
//  NcError error(NcError::verbose_nonfatal);
  NcError error(NcError::silent_nonfatal);
  try{

    //list of input files for which to calculate deviations
    std::string fileList;
    std::string fileOutList;
    //name of single input file
    std::string fileName;
    //name of single output file
    std::string outName;
    //file that holds averages
    std::string avgName;
    //Name of input variable
    std::string varName;
    //name of averaged input variable
    std::string avgVarName;
    std::string tname,latname,lonname,levname;
    std::string insuff,outsuff;
    bool is4D,ZtoGH,isHpa;
    std::string fourDval="F";
    std::string ZGHval="F";
    std::string Hpaval="F";
    std::string devName;
    std::string aDevName;
    std::string prevFile;
    std::string nextFile;
    bool PVCalc;
    bool GHCalc;
    size_t pos,len;
    std::string delim = ".";
    int tSteps;
    double missingNo;

    BeginCommandLine()
      CommandLineString(fileName, "in","");
      CommandLineString(outName, "out","");
      CommandLineString(fileList, "inlist", "");
      CommandLineString(fileOutList,"outlist","");
      CommandLineString(varName,"varname","");
      CommandLineString(avgName, "avg", "");
      CommandLineString(avgVarName, "avgname","");
      CommandLineBool(PVCalc,"pv");
      CommandLineBool(GHCalc,"z500");
      CommandLineBool(is4D,"is4D");
      CommandLineBool(ZtoGH,"gh");
      CommandLineBool(isHpa,"hpa");
      CommandLineString(tname,"tname","time");
      CommandLineString(levname,"levname","lev");
      CommandLineString(latname,"latname","lat");
      CommandLineString(lonname,"lonname","lon");
      CommandLineString(insuff,"insuff",".nc");
      CommandLineString(devName,"devname","");
      CommandLineString(aDevName,"adevname","");
      CommandLineString(outsuff,"outsuff","_devs.nc");
      CommandLineString(prevFile,"prev","");
      CommandLineString(nextFile,"next","");
      CommandLineInt(tSteps,"nt",0);
      CommandLineDouble(missingNo,"missingnum",1e20);
      ParseCommandLine(argc, argv);
    EndCommandLine(argv);
    AnnounceBanner();

    if ((!PVCalc) && (!GHCalc)){
      //Need to explicitly specify devName and ADevName
      if (devName==""){
        _EXCEPTIONT("Need to either specify PV calculation (--pv), Z500 calculation (--z500) or provide variable name (--devname).");
      }
      if (aDevName==""){
        _EXCEPTIONT("Need to either specify PV calculation (--pv), Z500 calculation (--z500) or provide variable name (--adevname).");
      }
    }
    if (fileName == "" && fileList == ""){
      _EXCEPTIONT("Need to specify either input file (--in) or file list (--inlist).");
    }
    if (fileName != "" && fileList != ""){
      _EXCEPTIONT("Cannot specify both file name (--in) and file list (--inlist).");
    }
    if (outName != "" && fileList != ""){
      _EXCEPTIONT("Currently cannot specify output name for list of files. This will only work with a single file (--in).");
    }

    if (ZtoGH){
      ZGHval="T";
    }
    if (PVCalc){
      devName = "DVPV";
			aDevName = "ADVPV";
    }else if (GHCalc){
      devName = "DZ";
			aDevName = "ADZ";
    }

   int nFiles,avgTime,nTime,nLat,nLon;
   double anomVal;
   if (PVCalc){
     anomVal = 1.3*std::pow(10,-6);
   }else if (GHCalc){
     anomVal = 170.;
   }
   //Create list of input files
    std::vector<std::string> InputFiles;
    std::vector<std::string> InputDevFiles;
    if (fileList != ""){
      GetInputFileList(fileList, InputFiles);
    }
    else{
      InputFiles.push_back(fileName);
    }
    nFiles = InputFiles.size();
    std::vector<std::string> OutputFiles;
  //Necessary because I need to write to a separate output directory!!
    if (fileOutList!=""){
      GetInputFileList(fileOutList,OutputFiles);
    }


    //Open averages file
    NcFile avgFile(avgName.c_str());
	if (!avgFile.is_valid()) {
		_EXCEPTION1("Cannot open NetCDF file \"%s\"", avgName.c_str());
	}

    //time data (average)
	NcDim *dimTime = avgFile.get_dim(tname.c_str());
	if (dimTime == NULL) {
		_EXCEPTION1("\"%s\" is missing dimension \"time\"", avgName.c_str());
	}
    avgTime = dimTime->size();
    NcVar *avgTimeVals = avgFile.get_var(tname.c_str());
   	if (avgTimeVals == NULL) {
		_EXCEPTION1("\"%s\" is missing variable \"time\"", avgName.c_str());
	}
    
    //averaged var
    NcVar *AvarData = avgFile.get_var(avgVarName.c_str());
   	if (AvarData == NULL) {
		_EXCEPTION2("\"%s\" is missing variable \"%s\"", avgName.c_str(), avgVarName.c_str());
	}
    int dim1 = AvarData->get_dim(0)->size();
    int dim2 = AvarData->get_dim(1)->size();
    int dim3 = AvarData->get_dim(2)->size();    
    //Open var files
    NcFile tempFile(InputFiles[0].c_str());
    NcDim *latDim = tempFile.get_dim(latname.c_str());
    NcDim *lonDim = tempFile.get_dim(lonname.c_str());
    nLat = latDim->size();
    nLon = lonDim->size();
    
    //Get 500 hPa level if 4D variable
    int pIndex = 10000000;
    if (is4D){
      fourDval="T";
      NcDim * lev = tempFile.get_dim(levname.c_str());
      int nLev = lev->size();
      NcVar *levvar = tempFile.get_var(levname.c_str());
      DataVector<double> pVec(nLev);
      levvar->get(&(pVec[0]),nLev);

      //Find the 500 mb level
      double pval = 50000.0;
      if (isHpa){
        Hpaval="T";
        pval = 500.0;
      }
      for (int x=0; x<nLev; x++){
        if (std::fabs(pVec[x]-pval)<0.0001){
          pIndex = x;
          break;
        }
      }
      if (pIndex > 999999){
        _EXCEPTIONT("Could not identify correct pressure level. Check file.");
      }
    }

    std::string strTimeUnits,strCalendar;
    if (tSteps<1){
      //Get info on number of time steps per day
      NcVar *tInfo = tempFile.get_var(tname.c_str());
      int tsz = tempFile.get_dim(tname.c_str())->size();
      DataVector<double> tvec(tsz);
      tInfo->set_cur((long)0);
      tInfo->get(&(tvec[0]),tsz);

      NcAtt *attCal = tInfo->get_att("calendar");
      if (attCal==NULL){
        strCalendar = "standard";
      }else{
        strCalendar = attCal->as_string(0);
      }
      NcAtt *attTime = tInfo->get_att("units");
      if (attTime == NULL){
        _EXCEPTIONT("Time variable has no units attribute.");
      }

      strTimeUnits = attTime->as_string(0);

      int tYear = 0;
      int tMonth = 0;
      int tDay = 0;
      int tHour = 0;
    
      ParseTimeDouble(strTimeUnits,strCalendar,tvec[0],tYear,\
        tMonth,tDay,tHour);

      int nextDay;
      for (int t=1; t<tsz; t++){
        ParseTimeDouble(strTimeUnits,strCalendar,tvec[t],tYear,\
          tMonth,nextDay,tHour);
        tSteps+=1;
        if (nextDay != tDay){
          break;
        }
      }
    }
    tempFile.close();
    std::cout<<"tsteps is "<<tSteps<<std::endl; 
		//First, the dev variable
    for (int x=0; x<nFiles; x++){

      NcFile infile(InputFiles[x].c_str());
      if(!infile.is_valid()){
        _EXCEPTION1("Unable to open file \"%s\".",InputFiles[0].c_str());
      }
      std::cout<<"Opening file "<<InputFiles[x].c_str()<<std::endl;
      NcDim *tDim = infile.get_dim(tname.c_str());
      nTime = tDim->size();

      NcVar *inTime = infile.get_var(tname.c_str());
      NcVar *inLat = infile.get_var(latname.c_str());
      NcVar *inLon = infile.get_var(lonname.c_str());
      NcVar *varData = infile.get_var(varName.c_str());

      

      int nDims = varData->num_dims();
      if (nDims > 3 && !is4D){
        _EXCEPTIONT("Error: variable has more than 3 dimensions, must use --is4D.");
      }

      NcAtt *attTime = inTime->get_att("units");
      if (attTime == NULL){
        _EXCEPTIONT("Time variable has no units attribute.");
      }

      strTimeUnits = attTime->as_string(0);
      NcAtt *attCal = inTime->get_att("calendar");
      if(attCal==NULL){
        strCalendar = "standard";
      }else{
        strCalendar = attCal->as_string(0);
      }


      DataVector<double> timeVals(nTime);
      inTime->set_cur((long) 0);
      inTime->get(&(timeVals[0]),nTime);

    //check if the file contains leap days 
      int leapYear = 0;
      int leapMonth = 0;
      int leapDay = 0;
      int leapHour = 0;

      DataVector <double> outputTime(nTime);

      int nLeapSteps = 0;
      int z=0;

      //How many leap day times step are there?
      if (strCalendar!="noleap" && strCalendar!="360_day"){
        for (int t=0; t<nTime; t++){
          ParseTimeDouble(strTimeUnits, strCalendar, timeVals[t], leapYear,\
            leapMonth, leapDay, leapHour);
          if ((leapMonth==2 && leapDay == 29)){
            nLeapSteps +=1;
          }     
          else{
            outputTime[z] = timeVals[t];
            z++;
          }
        }
      }else{
        outputTime = timeVals;
      }
	std::cout<<"DEBUG: THERE ARE "<<nLeapSteps<<" leap days in file"<<std::endl;
      
      //Create output file that corresponds to IPV data

      std::string strOutFile;
      if (fileOutList!=""){
        strOutFile = OutputFiles[x].c_str();
      }
      else if (outName != ""){
        strOutFile = outName;
      }else{
        pos = InputFiles[x].find(insuff);
        len = InputFiles[x].length();

        strOutFile = InputFiles[x].replace(pos,len,outsuff.c_str());
      }
      std::cout<<"Writing unsmoothed anomaly variables to file "<<strOutFile.c_str()<<std::endl;
			InputDevFiles.push_back(strOutFile);
      NcFile outfile(strOutFile.c_str(), NcFile::Replace, NULL,0,NcFile::Offset64Bits);
      if (!outfile.is_valid()){
	_EXCEPTION1("Couldn't open file %s",strOutFile.c_str());
      }
      int nOutTime;
      nOutTime = nTime-nLeapSteps;
      NcDim *tDimOut = outfile.add_dim(tname.c_str(),nOutTime);
      NcDim *latDimOut = outfile.add_dim(latname.c_str(),nLat);
      NcDim *lonDimOut = outfile.add_dim(lonname.c_str(),nLon);

      NcVar *tVarOut = outfile.add_var(tname.c_str(),ncDouble,tDimOut);
      CopyNcVarAttributes(inTime,tVarOut);      
      tVarOut->set_cur((long) 0);
      tVarOut->put(&(outputTime[0]),nOutTime);

      NcVar *latVarOut = outfile.add_var(latname.c_str(),ncDouble,latDimOut);
      copy_dim_var(inLat,latVarOut);
      NcVar *lonVarOut = outfile.add_var(lonname.c_str(),ncDouble,lonDimOut);
      copy_dim_var(inLon,lonVarOut);
      NcVar *devOut = outfile.add_var(devName.c_str(),ncDouble,tDimOut,latDimOut,lonDimOut);

      //Create variables for Deviations
      
      if (PVCalc){
	calcDevs(false,ZGHval,fourDval,pIndex, tSteps, nOutTime, strTimeUnits, strCalendar, \
        varData, devOut, AvarData, inTime, avgTimeVals, inLat,missingNo);
      }
      else if (GHCalc){
        calcDevs(true,ZGHval,fourDval,pIndex,tSteps, nOutTime, strTimeUnits, strCalendar,\
          varData, devOut,AvarData,inTime, avgTimeVals,inLat,missingNo);
      }
      else{
        calcDevs(false,ZGHval,fourDval,pIndex, tSteps, nOutTime, strTimeUnits, strCalendar, \
        varData, devOut, AvarData, inTime, avgTimeVals, inLat,missingNo);
      }

			outfile.close();
    }
		//Now, the smoothed devs
		int stepsLeft = 2*tSteps;
		int currMatIndex = 0;
	  double div = (double) 2*tSteps;
	  double invDiv = 1.0/div;
		//Create the two-day matrix which will hold the averaging window
		DataMatrix3D<double> twoDayMat(2*tSteps,nLat,nLon);
		DataMatrix3D<double> nextFileBuffer(tSteps,nLat,nLon);
		DataMatrix<double>aDevMat(nLat,nLon);
		int xcount=0;
                int prevYear,prevMonth,prevDay,prevHour,nextYear,nextMonth,nextDay,nextHour = 0;                
                bool isSequential;
		//Case 1: File contains less than two days' worth of data
		if (nTime<2*tSteps && stepsLeft>tSteps){
			while (stepsLeft>tSteps){
			    NcFile infile(InputDevFiles[xcount].c_str(),NcFile::Write);
			    if(!infile.is_valid()){
		   		   _EXCEPTION1("Unable to open file \"%s\".",InputDevFiles[xcount].c_str());
		 	   }
		  	  std::cout<<"Opening file "<<InputDevFiles[xcount].c_str()<<std::endl;
		   	  NcDim *tDim = infile.get_dim(tname.c_str());
		  	  nTime = tDim->size();
			    NcDim *latDim = infile.get_dim(latname.c_str());
			    NcDim *lonDim = infile.get_dim(lonname.c_str());

			  NcVar *varData = infile.get_var(devName.c_str());
		 	 NcVar *aDevOut = infile.add_var(aDevName.c_str(),ncDouble,tDim,latDim,lonDim);			
				
				varData->set_cur(0,0,0);
				aDevOut->set_cur(0,0,0);

				//Add the data to the two day matrix
				varData->get(&(twoDayMat[0][0][0]),nTime,nLat,nLon);
				aDevOut->put(&(twoDayMat[0][0][0]),nTime,nLat,nLon);
				stepsLeft -= nTime;
				currMatIndex+=nTime;
				xcount+=1;
			}		
		}
		//Matrix is either empty or at least filled with one day's worth of data
    for (int x=xcount; x<nFiles; x++){
			//Open the next file. Add a day's worth of data and start the moving window
			int tStart = 0;
      NcFile infile(InputDevFiles[x].c_str(),NcFile::Write);
      if(!infile.is_valid()){
        _EXCEPTION1("Unable to open file \"%s\".",InputDevFiles[x].c_str());
      }
      std::cout<<"Opening file "<<InputDevFiles[x].c_str()<<std::endl;
      NcDim *tDim = infile.get_dim(tname.c_str());
      NcVar *tVar = infile.get_var(tname.c_str());
      nTime = tDim->size();
	DataVector<double> tVarVals(nTime);
	tVar->set_cur((long) 0);
	tVar->get(&(tVarVals[0]),nTime);	
	NcDim *latDim = infile.get_dim(latname.c_str());
      NcDim *lonDim = infile.get_dim(lonname.c_str());
       NcAtt *attTime = tVar->get_att("units");
       if (attTime==NULL){
           _EXCEPTIONT("Time variable has no units attribute.");
       }
       strTimeUnits = attTime->as_string(0);
       NcAtt *attCal = tVar->get_att("calendar");
       if (attCal==NULL){
          strCalendar = "standard";
        }else{
          strCalendar = attCal->as_string(0);
        }
        if (strncmp(strCalendar.c_str(),"gregorian",9)==0){
          strCalendar = "standard";
        }
      ParseTimeDouble(strTimeUnits,strCalendar,tVarVals[nTime-1],prevYear,prevMonth,prevDay,prevHour);
      NcVar *varData = infile.get_var(devName.c_str());
      NcVar *aDevOut = infile.add_var(aDevName.c_str(),ncDouble,tDim,latDim,lonDim);
			
			//Load the buffer data from the next file (first day's worth of data)
                        //But first, check: are the files sequential?
			if (x<(nFiles-1)){
				NcFile nextFile(InputDevFiles[x+1].c_str());
				std::cout<<"DEBUG: next file is "<<InputDevFiles[x+1].c_str()<<std::endl;
                                NcVar * nextTvar = nextFile.get_var(tname.c_str());
                                int nextTsz = nextFile.get_dim(tname.c_str())->size();
				DataVector<double> nextVarVals(nextTsz);
				nextTvar->set_cur((long) 0);
				nextTvar->get(&(nextVarVals[0]),nextTsz);
				std::string strNextTime;
				std::string strNextCal;
				NcAtt * nextAttTime = nextTvar->get_att("units");
				if (nextAttTime==NULL){
					_EXCEPTIONT("Time variable has no units attribute.");
				}
				strNextTime = nextAttTime->as_string(0);

                                ParseTimeDouble(strNextTime,strCalendar,nextVarVals[0],nextYear,nextMonth,nextDay,nextHour);
                               std::cout<<"prev date is "<<prevYear<<"/"<<prevMonth<<"/"<<prevDay<<" and next date is "<<nextYear<<"/"<<nextMonth<<"/"<<nextDay<<std::endl; 
				isSequential = sequentialFiles(prevYear,prevMonth,prevDay,prevHour,\
                                  nextYear,nextMonth,nextDay,nextHour,strCalendar);
                                if (isSequential == true){
					NcVar *fillData = nextFile.get_var(devName.c_str());
					fillData->get(&(nextFileBuffer[0][0][0]),tSteps,nLat,nLon);
				}
                                else{
				std::cout<<"Missing some steps, filling in with buffer"<<std::endl;
				//Fill the buffer with fill values
					for (int t=0; t<tSteps; t++){
						for (int a=0; a<nLat; a++){
							for (int b=0; b<nLon; b++){
								nextFileBuffer[t][a][b] = -999999;
							}						
						}
					}					
				}
				nextFile.close();
			}
			
			varData->set_cur(0,0,0);
			aDevOut->set_cur(0,0,0);
			//Add the rest of the data
			if (stepsLeft>0){
				//If there is still data left to be filled, fill it starting at the mat index
				varData->get(&(twoDayMat[currMatIndex][0][0]),stepsLeft,nLat,nLon);
				//The starting calculation point will be at the one day point
				currMatIndex+=stepsLeft;
				tStart = stepsLeft-tSteps;
				//There should now be no steps left!
				stepsLeft-=stepsLeft;
			}
			if (tStart>0){
				//If t starts after 0, fill in the first few steps with the unsmoothed data which was just added
				aDevOut->put(&(twoDayMat[0][0][0]),tStart,nLat,nLon);
			}
			int buffIndex=0;
			for (int t=tStart; t<nTime;t++){
				aDevOut->set_cur(t,0,0);
				if (x==(nFiles-1) && t>=(nTime-tSteps)){
					varData->set_cur(t,0,0);
					varData->get(&(aDevMat[0][0]),1,nLat,nLon);
				}
				else{
					//Make sure that none of the two day mat contains fill data
					//if it does, the anomalies will be unsmoothed
					bool hasFill = false;
					for (int f=0; f<2*tSteps; f++){
						if (twoDayMat[f][0][0] <= -999998){
							hasFill = true;
							break;
						}
					}
					if (hasFill == true){
						varData->set_cur(t,0,0);
						varData->get(&(aDevMat[0][0]),1,nLat,nLon);
					}else{
						for (int a=0; a<nLat; a++){
							for (int b=0; b<nLon; b++){
								aDevMat[a][b] = 0.0;
								for (int n=0; n<2*tSteps; n++){
									aDevMat[a][b]+= (twoDayMat[n][a][b]*invDiv);
								}	
							}
						}
					}
					//Reset mat index if necessary
					if (currMatIndex>=2*tSteps){
						currMatIndex-=2*tSteps;
					}
					//Fill the two day mat with the next step 
					if (t<(nTime-tSteps)){
						varData->set_cur((t+tSteps),0,0);
						varData->get(&(twoDayMat[currMatIndex][0][0]),1,nLat,nLon);
					}
					else if(t>=(nTime-tSteps)){
						for (int a=0; a<nLat; a++){
							for (int b=0; b<nLon; b++){
								twoDayMat[currMatIndex][a][b] = nextFileBuffer[buffIndex][a][b];
							}
						}
						buffIndex+=1;
						if (buffIndex>=tSteps){
							buffIndex-=tSteps;
						}
					}	
				}
				aDevOut->put(&(aDevMat[0][0]),1,nLat,nLon);
				currMatIndex+=1;
			}	
	    std::cout<<"Finished writing smoothed anomalies to file "<<InputDevFiles[x].c_str()<<std::endl;
		}

  }
	catch (Exception & e){
    std::cout<<e.ToString() <<std::endl;
  }

}
