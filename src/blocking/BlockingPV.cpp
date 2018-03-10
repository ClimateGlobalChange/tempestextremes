/////////////////////////////////////////
///
///    \file blockingPV.cpp
//     \author Marielle Pinheiro
///    \version June 12, 2015
///


#include "BlockingUtilities.h"
#include "Interpolate.h"
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
#include <iostream>

//////////////////////////////////////////////////////////////////////// 
int main(int argc, char **argv){
  NcError error(NcError::verbose_nonfatal);

  try {

  //Input file (sometimes needs to be interpolated)
  std::string strfile_in;
  //Input file with PS
  std::string strfile_2d;
  //Input file list
  std::string fileList;
  //Input filelist with PS
  std::string fileList_PS;
  //axis names
  std::string tname,latname,lonname,levname;


  //Output file to be written
//  std::string strfile_out;
  //Interpolate, yes or no?
  bool interp_check;
  //If Pressure is in hPa, need to create new pressure axis
  bool is_hPa;
  //PV already calculated?
  bool has_PV;
  //If have PV,what is the name?
  std::string PVname, VPVname;
  //Name of U, V, T
  std::string uName, vName, tempName;
  //FIle name conventions
  std::string insuff,outsuff,insuff2d,outsuff2d;
  //Parse command line
  BeginCommandLine()

    CommandLineString(strfile_in, "in", "");
    CommandLineString(strfile_2d, "in2d", "");
    CommandLineString(fileList, "inlist","");
    CommandLineString(fileList_PS, "inlist2d","");
//    CommandLineString(strfile_out, "out", "");
    CommandLineBool(interp_check, "ipl");
    CommandLineBool(is_hPa, "hpa");
    CommandLineBool(has_PV, "hasPV");
    CommandLineString(PVname,"PVname","PV");
    CommandLineString(VPVname,"VPVname","VPV");
    CommandLineString(uName,"uname","U");
    CommandLineString(vName,"vname","V");
    CommandLineString(tempName,"tempname","T");
    CommandLineString(tname,"tname","time");
    CommandLineString(levname,"levname","lev");
    CommandLineString(latname,"latname","lat");
    CommandLineString(lonname,"lonname","lon");
    CommandLineString(insuff,"insuff",".nc");
    CommandLineString(outsuff,"outsuff","_integ.nc");
    CommandLineString(insuff2d,"insuff2d",".nc");
    CommandLineString(outsuff2d,"outsuff2d","_ipl.nc");
    ParseCommandLine(argc, argv);

  EndCommandLine(argv)
  AnnounceBanner();

  // Check input
/*  if (strfile_in == "") {
     _EXCEPTIONT("No input file (--in) specified");
  }*/
      if (strfile_in == "" && fileList == ""){
    _EXCEPTIONT("Need to specify either input file (--in) or file list (--inlist).");

    }
    if (strfile_in != "" && fileList != ""){
      _EXCEPTIONT("Cannot specify both file name (--in) and file list (--inlist).");
    }


/*
  if (strfile_out == "") {
     strfile_out = strfile_in;
     strfile_out = strfile_out.replace(strfile_out.end()-3,strfile_out.end(),"_integ.nc");
  }*/

  //Make list of files
  std::vector <std::string> InputFiles;
  if (fileList != ""){
    GetInputFileList(fileList,InputFiles);
  }else{
    InputFiles.push_back(strfile_in);
  }
  int nfiles = InputFiles.size();
  size_t pos, len;

  //If variables need interpolating, do this first!  
  if (interp_check) {

    if (strfile_2d == "" && fileList_PS == ""){
    _EXCEPTIONT("Need to specify either input file (--in2d) or file list (--inlist2d).");

    }
    if (strfile_2d != "" && fileList_PS != ""){
      _EXCEPTIONT("Cannot specify both file name (--in2d) and file list (--inlist2d).");
    }
    std::vector <std::string> PSFiles;
    if (fileList_PS != ""){
      GetInputFileList(fileList_PS,PSFiles);
    }else{
      PSFiles.push_back(strfile_2d);
    }
    
/*    if (strfile_2d == "") {
       _EXCEPTIONT("No input file (--in2d) specified");
    }*/

    int nfiles_PS = PSFiles.size();
    if (nfiles != nfiles_PS){
      _EXCEPTIONT("Number of files to be interpolated doesn't match 2D files.");
    }

    for (int x=0; x<nfiles_PS; x++){
      NcFile readin_int(InputFiles[x].c_str());
      if (!readin_int.is_valid()) {
        _EXCEPTION1("Unable to open file \"%s\" for reading",
          InputFiles[x].c_str());
      }
      //Check that time dimension is non-zero
      NcDim *time_int = readin_int.get_dim(tname.c_str());
      int time_len_int = time_int->size();
      if (time_len_int < 1){
        _EXCEPTIONT("This file has a time dimension length of 0.");
      }

    //Interpolated file name
      std::string strfile_copy = InputFiles[x];
      pos = strfile_copy.find(insuff2d);
      len = strfile_copy.length();
      std::string interp_file = strfile_copy.replace(pos,len,outsuff2d.c_str());
      std::cout<< "Interpolated file name is "<<interp_file<<std::endl;

    //open file that interpolated variables will be written to
      NcFile readin_out(interp_file.c_str(), NcFile::Replace, NULL,\
        0, NcFile::Offset64Bits);
      if (!readin_out.is_valid()) {
        _EXCEPTION1("Unable to open file \"%s\" for writing",
          interp_file.c_str());
      }

    //Interpolate variables and write to new file
      std::cout << "About to interpolate files. Entering interp util."<<std::endl;
//    interp_util(readin_int, strfile_2d, readin_out);
     
      std::string varlist;
      if (has_PV){
        varlist=PVname;
      }else{
        varlist=tempName+","+uName+","+vName;
      }
      interp_util(readin_int,strfile_2d,varlist,readin_out);
      readin_out.close();
      InputFiles[x] = interp_file;
    }
  }

  for (int x=0; x<nfiles; x++){
    NcFile readin(InputFiles[x].c_str());
    if (!readin.is_valid()) {
      _EXCEPTION1("Invalid file \"%s\"", InputFiles[x].c_str());
    }

  //Dimensions and associated variables
    NcDim *time = readin.get_dim(tname.c_str());
    int time_len = time->size();
    NcVar *timevar = readin.get_var(tname.c_str());

    NcDim *lev = readin.get_dim(levname.c_str());
    int lev_len = lev->size();
    NcVar *levvar = readin.get_var(levname.c_str());

    NcDim *lat = readin.get_dim(latname.c_str());
    int lat_len = lat->size();
    NcVar *latvar = readin.get_var(latname.c_str());

    NcDim *lon = readin.get_dim(lonname.c_str());
    int lon_len = lon->size();
    NcVar *lonvar = readin.get_var(lonname.c_str());
    
   //output file
    pos = InputFiles[x].find(insuff);
    len = InputFiles[x].length();

    //std::string strfile_out = "foo.nc"; 
    std::string strfile_out = InputFiles[x].replace(pos,len,outsuff.c_str());
//    std::string strfile_out = InputFiles[x].replace(InputFiles[x].end()-3,InputFiles[x].end(),"_integ.nc");
    NcFile file_out(strfile_out.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);

  //Dimensions: time, pressure, lat, lon
    NcDim *out_time = file_out.add_dim(tname.c_str(), time_len);
    NcDim *out_plev = file_out.add_dim(levname.c_str(), lev_len);
    NcDim *out_lat = file_out.add_dim(latname.c_str(), lat_len);
    NcDim *out_lon = file_out.add_dim(lonname.c_str(), lon_len);
 
  //COPY EXISTING DIMENSION VALUES
    NcVar *time_vals = file_out.add_var(tname.c_str(), ncDouble, out_time);
    NcVar *lev_vals = file_out.add_var(levname.c_str(), ncDouble, out_plev);
    NcVar *lat_vals = file_out.add_var(latname.c_str(), ncDouble, out_lat);
    NcVar *lon_vals = file_out.add_var(lonname.c_str(), ncDouble, out_lon); 


  //if Pressure axis is in hPa, need to change to Pa
    if (is_hPa){
      DataVector<double> levhPa(lev_len);
      levvar->set_cur((long) 0);
      levvar->get(&(levhPa[0]),lev_len);

      DataVector<double> levPa(lev_len);
      for (int p=0; p<lev_len; p++){
        levPa[p] = levhPa[p]*100.0;
      }
        lev_vals->set_cur((long) 0);
        lev_vals->put(&(levPa[0]),lev_len);
      }
      else if (!is_hPa){
        copy_dim_var(levvar, lev_vals);
      }

    copy_dim_var(timevar, time_vals);
    if (time_vals->get_att("calendar") == NULL){
      time_vals->add_att("calendar","standard");
    }
    copy_dim_var(latvar, lat_vals);
    copy_dim_var(lonvar, lon_vals);
    double lat_res,lon_res,dphi,dlambda,dp;;
    DataVector<double> coriolis(lat_len);
    DataVector<double> cosphi(lat_len);


    pv_vars_calc(lat_vals, lon_vals, lev_vals, lat_res, lon_res,\
      dphi, dlambda, dp, coriolis, cosphi);

    DataVector<double>pVec(lev_len);
    lev_vals->set_cur((long)0);
    lev_vals->get(&(pVec[0]),lev_len);

    DataMatrix3D<double> PVMat(lev_len,lat_len,lon_len);
    DataMatrix<double>IPVMat(lat_len,lon_len);
    NcVar * pv_var = NULL;
  //if PV hasn't already been calculated, calculate it!
   if (!has_PV){
    //Variables for calculations
      NcVar *temp = readin.get_var(tempName.c_str());
      if (temp == NULL){
        _EXCEPTION1("Could not find variable %s",tempName.c_str());
      }
      NcVar *uvar = readin.get_var(uName.c_str());
      if (uvar == NULL){
        _EXCEPTION1("Could not find variable %s",uName.c_str());
      } 
      NcVar *vvar = readin.get_var(vName.c_str());
      if (vvar == NULL){
        _EXCEPTION1("Could not find variable %s",vName.c_str());
      }

    //Add the PV variable to the file
      pv_var = file_out.add_var("PV", ncDouble, out_time, out_plev, out_lat, out_lon);
    


      DataMatrix3D<double> PTVar(lev_len,lat_len,lon_len);
      DataMatrix3D<double> TVar(lev_len,lat_len,lon_len);
      DataMatrix3D<double>RVVar(lev_len,lat_len,lon_len);
      DataMatrix3D<double>UMat(lev_len,lat_len,lon_len);
      DataMatrix3D<double>VMat(lev_len,lat_len,lon_len);
   //CALCULATE PV AND IPV PER TIME STEP
      for (int t=0; t<time_len; t++){
      //Calculate PT
        temp->set_cur(t,0,0,0);
        temp->get(&(TVar[0][0][0]),1,lev_len,lat_len,lon_len);
        PT_calc(lev_len,lat_len,lon_len,TVar, lev_vals, PTVar);

      //Calculate relative vorticity 

        uvar->set_cur(t,0,0,0);
        uvar->get(&(UMat[0][0][0]),1,lev_len,lat_len,lon_len);

        vvar->set_cur(t,0,0,0);
        vvar->get(&(VMat[0][0][0]),1,lev_len,lat_len,lon_len);
        rVort_calc(lev_len,lat_len,lon_len,UMat,VMat, dphi, dlambda, cosphi, RVVar);

        PV_calc(lev_len,lat_len,lon_len, UMat, VMat, PTVar, RVVar, pVec, coriolis,cosphi, dphi, dlambda,\
          lat_res, lon_res, PVMat);

        pv_var->set_cur(t,0,0,0);
        pv_var->put(&(PVMat[0][0][0]),1,lev_len,lat_len,lon_len);
      }
    }
    if (has_PV){
      pv_var = readin.get_var(PVname.c_str());
    }
    
    NcVar *intpv_var = file_out.add_var(VPVname.c_str(), ncDouble, out_time, out_lat, out_lon);
    for (int t=0; t<time_len; t++){
      pv_var->set_cur(t,0,0,0);
      pv_var->get(&(PVMat[0][0][0]),1,lev_len,lat_len,lon_len);
      IPV_calc(lev_len,lat_len,lon_len,lat_res,pVec,PVMat,IPVMat);
      intpv_var->set_cur(t,0,0);
      intpv_var->put(&(IPVMat[0][0]),1,lat_len,lon_len);

    }
/*
    NcVar *avgt_var = file_out.add_var("AVGT", ncDouble, out_time, out_lat, out_lon);
    NcVar *avgu_var = file_out.add_var("AVGU", ncDouble, out_time, out_lat, out_lon);
    NcVar *avgv_var = file_out.add_var("AVGV", ncDouble, out_time, out_lat, out_lon);
  
    VarPressureAvg(temp,lev_vals,avgt_var);
    VarPressureAvg(uvar,lev_vals,avgu_var);
    VarPressureAvg(vvar,lev_vals,avgv_var);
*/


    std::cout<<"About to close files."<<std::endl;
 //Close input files
    readin.close();

  //Close output file
    file_out.close();
  }
  } catch (Exception & e) {
    std::cout << e.ToString() << std::endl;
  }
}
