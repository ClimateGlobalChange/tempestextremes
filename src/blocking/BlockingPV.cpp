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
  //Output file to be written
  std::string strfile_out;
  //Interpolate, yes or no?
  bool interp_check;
  //If Pressure is in hPa, need to create new pressure axis
  bool is_hPa;
  
  //Parse command line
  BeginCommandLine()

    CommandLineString(strfile_in, "in", "");
    CommandLineString(strfile_2d, "in2d", "");
    CommandLineString(strfile_out, "out", "");
    CommandLineBool(interp_check, "ipl");
    CommandLineBool(is_hPa, "hpa");
    ParseCommandLine(argc, argv);

  EndCommandLine(argv)
  AnnounceBanner();

  // Check input
  if (strfile_in == "") {
     _EXCEPTIONT("No input file (--in) specified");
  }
  if (strfile_out == "") {
     strfile_out = strfile_in;
     strfile_out = strfile_out.replace(strfile_out.end()-3,strfile_out.end(),"_integ.nc");
  }

  //If variables need interpolating, do this first!  
  if (interp_check) {
    if (strfile_2d == "") {
       _EXCEPTIONT("No input file (--in2d) specified");
    }
    NcFile readin_int(strfile_in.c_str());
    if (!readin_int.is_valid()) {
      _EXCEPTION1("Unable to open file \"%s\" for reading",
        strfile_in.c_str());
    }
    //Check that time dimension is non-zero
    NcDim *time_int = readin_int.get_dim("time");
    int time_len_int = time_int->size();
    if (time_len_int < 1){
      _EXCEPTIONT("This file has a time dimension length of 0.");
    }

    //Interpolated file name
    std::string interp_file = strfile_in.replace(strfile_in.end()-3,strfile_in.end(),"_ipl.nc");
    std::cout<< "Interpolated file name is "<<interp_file<<std::endl;
    std::cout<<"About to open interpolated file."<<std::endl;
    //open interpolated file
    NcFile readin_out(interp_file.c_str(), NcFile::Replace, NULL,\
      0, NcFile::Offset64Bits);
    if (!readin_out.is_valid()) {
      _EXCEPTION1("Unable to open file \"%s\" for writing",
        interp_file.c_str());
    }

    //Interpolate variable and write to new file
    std::cout << "About to interpolate files. Entering interp util."<<std::endl;
    interp_util(readin_int, strfile_2d, readin_out);
    readin_out.close();
    strfile_in = interp_file;
  }

  NcFile readin(strfile_in.c_str());
  if (!readin.is_valid()) {
    _EXCEPTION1("Invalid file \"%s\"", strfile_in.c_str());
  }

  //Dimensions and associated variables
  NcDim *time = readin.get_dim("time");
  int time_len = time->size();
  NcVar *timevar = readin.get_var("time");

  NcDim *lev = readin.get_dim("lev");
  int lev_len = lev->size();
  NcVar *levvar = readin.get_var("lev");

  NcDim *lat = readin.get_dim("lat");
  int lat_len = lat->size();
  NcVar *latvar = readin.get_var("lat");

  NcDim *lon = readin.get_dim("lon");
  int lon_len = lon->size();
  NcVar *lonvar = readin.get_var("lon");
    
  //Variables for calculations
  NcVar *temp = readin.get_var("T");
  NcVar *uvar = readin.get_var("U");
  NcVar *vvar = readin.get_var("V");

  //Create output file
  NcFile file_out(strfile_out.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);

  //Dimensions: time, pressure, lat, lon
  NcDim *out_time = file_out.add_dim("time", time_len);
  NcDim *out_plev = file_out.add_dim("lev", lev_len);
  NcDim *out_lat = file_out.add_dim("lat", lat_len);
  NcDim *out_lon = file_out.add_dim("lon", lon_len);
 
  //COPY EXISTING DIMENSION VALUES
  NcVar *time_vals = file_out.add_var("time", ncDouble, out_time);
  NcVar *lev_vals = file_out.add_var("lev", ncDouble, out_plev);
  NcVar *lat_vals = file_out.add_var("lat", ncDouble, out_lat);
  NcVar *lon_vals = file_out.add_var("lon", ncDouble, out_lon); 


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

  //Data for PT, RV, PV calculations
  DataVector<double> coriolis(lat_len);
  DataVector<double> cosphi(lat_len);

  double lat_res;
  double lon_res;

  double dphi;
  double dlambda;
  double dp;

  pv_vars_calc(lat_vals, lon_vals, lev_vals, lat_res, lon_res,\
    dphi, dlambda, dp, coriolis, cosphi);

//Calculate PT
  DataMatrix4D<double> PTVar(time_len,lev_len,lat_len,lon_len);
  PT_calc(temp, lev_vals, PTVar);

  //Calculate relative vorticity and add to outfile
  DataMatrix4D<double>RVVar(time_len,lev_len,lat_len,lon_len);
  rVort_calc(uvar, vvar, dphi, dlambda, cosphi, RVVar);

  NcVar *pv_var = file_out.add_var("PV", ncDouble, out_time, out_plev, out_lat, out_lon);
  NcVar *intpv_var = file_out.add_var("IPV", ncDouble, out_time, out_lat, out_lon);
  PV_calc(uvar, vvar, PTVar, RVVar, lev_vals, coriolis,cosphi, dphi, dlambda,\
    lat_res, lon_res, pv_var, intpv_var);

  std::cout<<"About to close files."<<std::endl;
 //Close input files
  readin.close();

  //Close output file
  file_out.close();

  } catch (Exception & e) {
    std::cout << e.ToString() << std::endl;
  }
}
