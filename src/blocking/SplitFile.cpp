//////////////////////////
///
///      \file split_file.cpp
///      \author Marielle Pinheiro
///      \version May 11,2016


#include "BlockingUtilities.h"
#include "NetCDFUtilities.h"
#include "Exception.h"
#include "CommandLine.h"
#include "netcdfcpp.h"
#include "DataArray1D.h"
#include "DataArray3D.h"
#include <cstdlib>
#include <cmath>
#include "TimeObj.h"
#include <string>
#include <vector>

int vspfunc(char * buffer, char *format, ...)
{
   va_list aptr;
   int ret;

   va_start(aptr, format);
   ret = vsprintf(buffer, format, aptr);
   va_end(aptr);

   return(ret);
}

void replace_date(std::string str,
                         int y1,int y2,
                         int m1,int m2,
                         int d1,int d2,
                         std::string &f1){

  char fmt[] = "%04d-%02d-%02d-00000";
  char to_replace[24];
  char replace_with[24];
  
  vspfunc(to_replace, fmt, y1, m1, d1);
  vspfunc(replace_with, fmt, y2, m2, d2);

  std::string str1(to_replace);
  std::string str2(replace_with);

  size_t f = str.find(str1);
  f1=str.replace(f,str1.length(),str2);
}

void add_vars_to_file(NcFile & infile,
                      NcFile & outfile,
                      int ndims,
                      int newTLen1,
                      int pos0,
                      std::string timename,
                      std::string levname,
                      std::string latname,
                      std::string lonname,
                      std::vector<std::string> vars){
  NcVar * timeVar = infile.get_var(timename.c_str());
  NcVar * latVar = infile.get_var(latname.c_str());
  int nLat = latVar ->get_dim(0)->size();
  NcVar * lonVar = infile.get_var(lonname.c_str());
  NcVar * levvar = infile.get_var(levname.c_str());
  NcVar * lev1var = NULL;
  NcDim * lev1 = NULL;
  int nLon = lonVar ->get_dim(0)->size();
  int nLev = 0;
  if (levvar != NULL){
    nLev = levvar->get_dim(0)->size();
    std::cout<<"nLev is "<<nLev<<std::endl;
  }

  DataArray1D<double> t1(newTLen1);
  timeVar->set_cur((long) pos0);
  timeVar->get(&(t1[0]),newTLen1);


  NcDim * time1 = outfile.add_dim(timename.c_str(), newTLen1);
  if (nLev >0){
    NcDim * lev1 = outfile.add_dim(levname.c_str(),nLev);
  }
  NcDim * lat1 = outfile.add_dim(latname.c_str(), nLat);
  NcDim * lon1 = outfile.add_dim(lonname.c_str(),nLon);

  NcVar * t1var = outfile.add_var(timename.c_str(), ncDouble, time1);
  t1var->set_cur((long) 0);
  t1var->put(&(t1[0]),newTLen1);
  CopyNcVarAttributes(timeVar,t1var);

  if (nLev > 0){
    NcVar * lev1var = outfile.add_var(levname.c_str(), ncDouble,lev1);
    copy_dim_var(levvar,lev1var);
  }

  NcVar * lat1var = outfile.add_var(latname.c_str(), ncDouble,lat1);
  copy_dim_var(latVar,lat1var);
  NcVar * lon1var = outfile.add_var(lonname.c_str(),ncDouble,lon1);
  copy_dim_var(lonVar,lon1var);

  if (nLev>0){
    DataArray4D<double> vMat(newTLen1,nLev,nLat,nLon);
    for (int v=0; v<vars.size(); v++){
      NcVar * invar = infile.get_var(vars[v].c_str());
      invar->set_cur(pos0,0,0,0);
      std::cout<<"Setting t to "<<pos0<<std::endl;
      invar->get(&(vMat[0][0][0][0]),newTLen1,nLev,nLat,nLon);
      NcVar * outvar = outfile.add_var(vars[v].c_str(), ncDouble, time1,lev1,lat1,lon1);
      outvar->set_cur(0,0,0,0);
      outvar->put(&(vMat[0][0][0][0]),newTLen1,nLev,nLat,nLon);
    }
  }else{
    DataArray3D<double> vMat(newTLen1,nLat,nLon);
    for (int v=0; v<vars.size(); v++){
      NcVar * invar = infile.get_var(vars[v].c_str());
      invar->set_cur(pos0,0,0);
      invar->get(&(vMat[0][0][0]),newTLen1,nLat,nLon);
      NcVar * outvar = outfile.add_var(vars[v].c_str(), ncDouble, time1,lat1,lon1);
      outvar->set_cur(0,0,0);
      outvar->put(&(vMat[0][0][0]),newTLen1,nLat,nLon);
    }
  }
}



int main(int argc, char ** argv ){
  NcError error(NcError::verbose_nonfatal);
  try{
  std::string fName;
  std::string fOut1;
  std::string fOut2;
  bool rename_f;
  std::string varlist;
  int splitYear, splitMonth, splitDay, splitHour;
  std::string timename,levname,latname,lonname;

  BeginCommandLine()
    CommandLineString(fName,"in","");
    CommandLineString(fOut1,"out1","");
    CommandLineString(fOut2,"out2","");
//    CommandLineBool(rename_f,"rename");
    CommandLineString(varlist,"vars","");
    CommandLineInt(splitYear,"year",0);
    CommandLineInt(splitMonth,"month",0);
    CommandLineInt(splitDay,"day",0);
    CommandLineInt(splitHour,"hour",0);
    CommandLineString(timename,"timename","time");
    CommandLineString(levname,"levname","lev");
    CommandLineString(latname,"latname","lat");
    CommandLineString(lonname,"lonname","lon");

    ParseCommandLine(argc,argv);
  EndCommandLine(argv)

  if (varlist == ""){
    _EXCEPTIONT("Need to provide list of variables to split (--vars).");
  }

/*  if ((!rename_f) && (fOut1 == "")){
    _EXCEPTIONT("Need to either specify --rename or provide a new filename for out1.");
  }*/
  if (fOut1 == "" && fOut2 == ""){
    _EXCEPTIONT("Need to provide output file names for the two new files (--out1 and --out2).");
  }

  //Split var list  
  std::string delim = ",";
  size_t pos = 0;
  std::string token;
  std::vector<std::string> varVec;

  while((pos = varlist.find(delim)) != std::string::npos){
    token = varlist.substr(0,pos);
    varVec.push_back(token);
    varlist.erase(0,pos + delim.length());
  }
  varVec.push_back(varlist);
  for (int v=0; v<varVec.size(); v++){
    std::cout<<"Vector contains string "<<varVec[v].c_str()<<std::endl;
  }

 

  //Open the input file and the time variable
  NcFile infile(fName.c_str());
  int ndims = infile.num_dims();
  std::cout<<"Input file has "<<ndims<<" dimensions."<<std::endl;

  NcVar *timeVar = infile.get_var(timename.c_str());
  NcVar *latVar = infile.get_var(latname.c_str());
  NcVar *lonVar = infile.get_var(lonname.c_str());
  NcVar *levVar = infile.get_var(levname.c_str());  


  int nTime = infile.get_dim(timename.c_str())->size();
  int nLat = infile.get_dim(latname.c_str())->size();
  int nLon = infile.get_dim(lonname.c_str())->size();
  int nLev = 0;
  if (levVar != NULL){
    nLev = infile.get_dim(levname.c_str())->size();
  }

  DataArray1D<double> timeVec(nTime);
  timeVar ->set_cur((long) 0);
  timeVar ->get(&(timeVec[0]),nTime);

 
  //Find the index that has a date with a day 1
  NcAtt *attTime = timeVar->get_att("units");
  std::string strTimeUnits = attTime->as_string(0);
  NcAtt *attCal = timeVar->get_att("calendar");
  std::string strCalendar = attCal ->as_string(0);

  int dateYear = 0;
  int dateMonth = 0;
  int dateDay = 0;
  int dateHour = 0;

  int cutIndex = -1000;

  for (int t=0; t<nTime; t++){
    ParseTimeDouble(strTimeUnits,strCalendar,timeVec[t],\
      dateYear,dateMonth,dateDay,dateHour);
  /*  if (dateDay == 1 && dateHour == 0){
      monthStartIndex = t;
      std::cout<< "found month start index at t="<<t<<std::endl;
      break;
    }*/
    if (dateYear == splitYear && dateMonth == splitMonth && \
          dateDay == splitDay && dateHour == splitHour){
      cutIndex = t;
      break;
    }
  }
  std::cout<<"Cut index is "<<cutIndex<<std::endl;
  if (cutIndex<0){
    _EXCEPTIONT("Check file-- no beginning month date found.");
  }  
  //use ncks to split the files
  std::string cmd1 = "ncks -d " + timename + ",0," + std::to_string((long long)(cutIndex-1)) + " " + fName + " " + fOut1;
  std::string cmd2 = "ncks -d " + timename + ","+ std::to_string((long long)(cutIndex)) + ", " + fName + " " + fOut2;
  system(cmd1.c_str());
  system(cmd2.c_str());

// std::cout<<cmd1.c_str()<<std::endl;
// std::cout<<cmd2<<std::endl;
  //create the output file names
 /* int year1=0;
  int month1=0;
  int day1=0;
  int hour1=0;

  int year2=0;
  int month2=0;
  int day2=1;
  int hour2=0;

  ParseTimeDouble(strTimeUnits,strCalendar,timeVec[0],\
    year1,month1,day1,hour1);

  if (month1==12){
    year2=year1+1;
    month2=1;
  }else{
    year2=year1;
    month2=month1+1;
  }
*/


//Deal with the file names for the new and old files
  
/*  if (rename_f){
    std::string fNew;

    std::string fname_copy = fName;
    fNew = fname_copy.replace(fname_copy.end()-3,fname_copy.end(), "_1.nc");
   //execute a system command to change the old file name
    std::string cmd = "mv " + fName + " " + fNew;
    //std::cout << cmd.c_str()<<std::endl;
    system(cmd.c_str());
  }

  if (fOut2 == ""){
    replace_date(fName,year1,year2,month1,month2,day1,day2,fOut2);
  }
  if (fOut1 == ""){
    fOut1 = fName;
  }*/

  //Dimensions of split time
/*  int newTLen1 = cutIndex;
  int newTLen2 = nTime-newTLen1;
  std::cout<<"New time dimensions are "<<newTLen1<<" and "<<newTLen2<<std::endl;

  //new time vectors
  DataArray1D<double> t1(newTLen1);
  timeVar->set_cur((long) 0);
  timeVar->get(&(t1[0]),newTLen1);

  DataArray1D<double> t2(newTLen2);
  timeVar->set_cur((long) cutIndex);
  timeVar->get(&(t2[0]),newTLen2);


  //Now create the new files
  NcFile file1(fOut1.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);
  add_vars_to_file(infile,file1, ndims,newTLen1,0,timename,levname,latname,lonname,varVec);
  NcFile file2(fOut2.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);
  add_vars_to_file(infile,file2,ndims,newTLen2,newTLen1,timename,levname,latname,lonname,varVec);
*/

  }
  catch(Exception &e){
    std::cout<<e.ToString()<<std::endl;
  }
}
