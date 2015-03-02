/////////////////////////////////////////
///
///    \file CLIVAR_blocks.cpp
//     \author Marielle Pinheiro
///    \version February 3, 2015
///


#include "CLIVAR_block_utilities.h"
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

  //Input file with variable to be interpolated
  std::string strfile_in;
  //Input file with PS
  std::string strfile_2d;
  //Output file to be written
  std::string strfile_out;
  
  //Parse command line
  BeginCommandLine()

    CommandLineString(strfile_in, "in", "");
    CommandLineString(strfile_2d, "in2d", "");
    CommandLineString(strfile_out, "out", "");

    ParseCommandLine(argc, argv);

  EndCommandLine(argv)
  AnnounceBanner();

    // Check input
    if (strfile_in == "") {
       _EXCEPTIONT("No input file (--in) specified");
    }
    if (strfile_2d == "") {
       _EXCEPTIONT("No input file (--in2d) specified");
    }
    if (strfile_out == "") {
       _EXCEPTIONT("No output file (--out) specified");
    }

  //Open input file (as read-only)
  NcFile readin(strfile_in.c_str());
  std::cout<< "Reading file "<< strfile_in.c_str() << std::endl;  

  //Open 2D input file
  NcFile readin_2d(strfile_2d.c_str());
  std::cout<< "Reading file " << strfile_2d.c_str() << std::endl;
  
  //Dimensions and associated variables
  NcDim *time = readin.get_dim("time");
  int time_len = time->size();
  NcVar *timevar = readin.get_var("time");

  NcDim *lev = readin.get_dim("lev");
  int lev_len = lev->size();
  
  NcDim *lat = readin.get_dim("lat");
  int lat_len = lat->size();
  NcVar *latvar = readin.get_var("lat");

  NcDim *lon = readin.get_dim("lon");
  int lon_len = lon->size();
  NcVar *lonvar = readin.get_var("lon");


  //Coefficients
  NcVar *hyam = readin.get_var("hyam");
  NcVar *hybm = readin.get_var("hybm");

  //Variables
  NcVar *temp = readin.get_var("T");
  NcVar *uvar = readin.get_var("U");
  NcVar *vvar = readin.get_var("V");
  NcVar *ps = readin_2d.get_var("PS");

  //Create pressure level vector
  int plev_len = (100000.0-5000.0)/(5000.0);

  DataVector<double> pVals(plev_len);

  for (int i=0; i<plev_len; i++){
    double pNum = 5000.0 * (i+1);
    pVals[i] = pNum;
  }

  //Create output file
  NcFile file_out(strfile_out.c_str(), NcFile::Replace);

  //Dimensions: time, pressure, lat, lon
  NcDim *out_time = file_out.add_dim("time", time_len);
  NcDim *out_plev = file_out.add_dim("lev", plev_len);
  NcDim *out_lat = file_out.add_dim("lat", lat_len);
  NcDim *out_lon = file_out.add_dim("lon", lon_len);
 
  //COPY EXISTING DIMENSION VALUES
  NcVar *time_vals = file_out.add_var("time", ncDouble, out_time);
  NcVar *lat_vals = file_out.add_var("lat", ncDouble, out_lat);
  NcVar *lon_vals = file_out.add_var("lon", ncDouble, out_lon); 

  copy_dim_var(timevar, time_vals);
  copy_dim_var(latvar, lat_vals);
  copy_dim_var(lonvar, lon_vals);

  //Give pressure level dimension values
  NcVar *plev_vals = file_out.add_var("lev", ncDouble, out_plev);
  plev_vals-> set_cur((long) 0);
  plev_vals->put(&(pVals[0]), plev_len);

  //Add interpolated variables to outfile
  NcVar *out_temp = file_out.add_var("T", ncDouble, out_time, out_plev, out_lat, out_lon);
  interpolate_lev(temp, hyam, hybm, ps, plev_vals, out_temp);

  NcVar *out_uvar = file_out.add_var("U", ncDouble, out_time, out_plev, out_lat, out_lon);
  interpolate_lev(uvar, hyam, hybm, ps, plev_vals, out_uvar);

  NcVar *out_vvar = file_out.add_var("V", ncDouble, out_time, out_plev, out_lat, out_lon);
  interpolate_lev(vvar, hyam, hybm, ps, plev_vals, out_vvar);

  ////////////////////////////////////////////////////////////

  //Variables for PV calculation
//  double radius = 6371000.0;
//  double pi = 4.0*cmath::atan(1.0);
//  double radian = 180.0/pi;
//  double sigma = cmath::pow(7.2921, -5.0);

  DataVector<double> coriolis(lat_len);
  DataVector<double> cosphi(lat_len);

  double dphi;
  double dlambda;
  double dp;

  pv_vars_calc(lat_vals, lon_vals, plev_vals, dphi, dlambda, dp, coriolis, cosphi);

  //Calculate PT and add to outfile
  NcVar *pt_var = file_out.add_var("PT", ncDouble, out_time, out_plev, out_lat, out_lon);
  PT_calc(out_temp, plev_vals, pt_var);

  //Calculate relative vorticity and add to outfile
  NcVar *rvort_var = file_out.add_var("REL_VORT", ncDouble, out_time, out_plev, out_lat, out_lon);
  rVort_calc(out_uvar, out_vvar, dphi, dlambda, cosphi, rvort_var);
 

 //Close input files
  readin.close();
  readin_2d.close();

  //Close output file
  file_out.close();
}
