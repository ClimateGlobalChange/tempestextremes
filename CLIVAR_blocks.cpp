/////////////////////////////////////////
///
///    \file CLIVAR_blocks.cpp
//     \author Marielle Pinheiro
///    \version February 3, 2015
///


#include "CLIVAR_block_utilities.h"
#include "interpolate.h"
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

  //Input file (sometimes needs to be interpolated)
  std::string strfile_in;
  //Input file with PS
  std::string strfile_2d;
  //Output file to be written
  std::string strfile_out;
  //Interpolate, yes or no?
  std::string interp_check;
  
  //Parse command line
  BeginCommandLine()

    CommandLineString(strfile_in, "in", "");
    CommandLineString(strfile_2d, "in2d", "");
    CommandLineString(strfile_out, "out", "");
    CommandLineString(interp_check, "ipl", "");
    ParseCommandLine(argc, argv);

  EndCommandLine(argv)
  AnnounceBanner();

    // Check input
    if (strfile_in == "") {
       _EXCEPTIONT("No input file (--in) specified");
    }
    if (strfile_out == "") {
       _EXCEPTIONT("No output file (--out) specified");
    }

  //If variables need interpolating, do this first!i  
  if (interp_check == "y" or interp_check == "yes"){
    if (strfile_2d == "") {
       _EXCEPTIONT("No input file (--in2d) specified");
    }
    //Interpolated file name
    std::string interp_file = "test_interp.nc";

    //Interpolate variable and write to new file
    NcFile readin_int(strfile_in.c_str());
    interp_util(readin_int, strfile_2d, interp_file);
   // std::cout<<"About to close file."<<std::endl;
   // readin_int.close();
    //Replace infile name with interpolated file name
    strfile_in = interp_file;
    std::cout<<"Input file is now "<<strfile_in<<std::endl;
  }
  //Read input file
  NcFile readin(strfile_in.c_str());

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

  copy_dim_var(timevar, time_vals);
  copy_dim_var(levvar, lev_vals);
  copy_dim_var(latvar, lat_vals);
  copy_dim_var(lonvar, lon_vals);

  //Data for PT, RV, PV calculations
  DataVector<double> coriolis(lat_len);
  DataVector<double> cosphi(lat_len);

  double dphi;
  double dlambda;
  double dp;

  pv_vars_calc(lat_vals, lon_vals, lev_vals, dphi, dlambda, dp, coriolis, cosphi);
  std::cout<<"dphi: "<<dphi<<" dlambda: "<<dlambda<<" dp: "<<dp<<std::endl;

  //Calculate PT and add to outfile
  NcVar *pt_var = file_out.add_var("PT", ncDouble, out_time, out_plev, out_lat, out_lon);
  PT_calc(temp, lev_vals, pt_var);

  //Calculate relative vorticity and add to outfile
  NcVar *rvort_var = file_out.add_var("REL_VORT", ncDouble, out_time, out_plev, out_lat, out_lon);
  rVort_calc(uvar, vvar, dphi, dlambda, cosphi, rvort_var);

  NcVar *pv_var = file_out.add_var("PV", ncDouble, out_time, out_plev, out_lat, out_lon);
  NcVar *intpv_var = file_out.add_var("IPV", ncDouble, out_time, out_lat, out_lon);
  PV_calc(uvar, vvar, pt_var, rvort_var, lev_vals, coriolis,cosphi, dphi, dlambda, pv_var, intpv_var);

 //Close input files
  readin.close();

  //Close output file
  file_out.close();
}
