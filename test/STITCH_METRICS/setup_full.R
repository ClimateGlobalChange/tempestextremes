#This is an example namelist for an output report
#The requirements are as follows:

#1) StitchBlobs files from each dataset
#2) Corresponding BlobStats files from each dataset
#3) Optional but recommended: BlobStats output from
# the separate DetectBlobs function
# (This will help deal with cases where two blobs
# merge into a single blob at a later time)

#REPORT HEADER: ANY TEXT THAT YOU WANT TO GO AT THE TOP OF THE REPORT
#NOTE: must be in quotes. Careful of any special characters
metadata_datasets="This is the results summary for DJF NP (June 1980-August 2004). 
There are 4 reanalysis datasets (ERA-Interim, MERRA2, JRA, and CFSR) and 12 historical model outputs. 
All of the reanalysis data is 6-hourly, and all of the model data is daily (12Z). \n\n
Intercomparison analysis for probability of co-occurrence, spatial similarity, Pearson pattern correlation, 
and RMSE is performed using only the common subset of timesteps, and all data has
been regridded to 1x1 degree (will assess whether changing the regrid resolution affects 
intercomparison results).
The individual dataset information for number of blobs, number of blocked days,
etc retains all time steps."
#This is the number of years in your dataset! (Integer)
nyears<-2004-1980+1

###########
#FILE INFO
#Will always output RData files
#But there is an option to also output text files
# or CSV files with the table data
#Set to TRUE if you wish to have one or both of these
output_txt<-FALSE
output_csv<-FALSE

#Use DetectBlobs inputs? (TRUE/FALSE)
use_detectblob<-TRUE

#Names of the datasets, which will be used to differentiate between them
#in the tables (Vector of strings)
#Varnames<-c("MERRA","JRA")
Varnames<-c("CMCC-CESM","CMCC-CM","CNRM-CM5","CanESM2","GFDL-CM3","GFDL-ESM2M",
           "MIROC-ESM","MIROC5","MPI-ESM-MR","MRI-CGCM3",
           "MRI-ESM1","NorESM1-M","CFSR","ERA","MERRA","JRA")

#The spatial resolutions of the datasets (Vector of strings. Optional-- delete this variable if 
# you don't want to use it)
resolutions<-c("3.3341x3.75","TM159 (0.75)","T85 (1.39)","T42 (2.79)",
               "2x2.5","2x2.5","T42 (2.8)","T85 (1.39)","T63 (1.88)",
               "T106 (1.12)","T106 (1.12)","1.9x2.5","0.5x0.5","1x1",
               "0.5x0.625","1.25x1.25")
#Working directory (where all of the R functions are)
work_dir<-"~/tempestextremes/test/STITCH_METRICS"
#Output directory (where all files will go)
output_dir<-"~/BLOBSTATS_FILES"
#Name of output subdirectory
#This will be created in the output directory
output_subdir<-"DJF_NP"
#Name of the file prefix. Output file will be [prefix]_[Varname]_[suffix] 
output_prefix<-"DJF_NP"

#Output name for the master namelist that will be generated
fname_namelist<-"DJF_NP_namelist_master_test.R"

#NOTE THAT ALL OF THE FILES CONTAINED IN THESE LISTS SHOULD HAVE FULL PATHNAMES
#Note that the input directory is appended to these list files!


#List of BlobStats files (from StitchBlobs output) to read into readfiles.R
#USER_DEFINED
#Make the lists on the fly using the search strings above!
##########################

#USER_DEFINED
#Search string for stitch files
stitch_searchstring<-"*DJF_NP*stats.txt"
nostitch_searchstring<-"*DJF_NP*stats_nostitch.txt"
netcdf_searchstring<-"*DJF_NP*blobs.nc"
netcdf2_searchstring<-"*DJF*blobs.nc"


#Input directory
input_directory<-"~/BLOBSTATS_FILES"
stitchlist_fnames<-sprintf("%s/%s/%s_stitch_list_NA",input_directory,Varnames,Varnames)
nostitch_fnames<-sprintf("%s/%s/%s_nostitch_list_NA",input_directory,Varnames,Varnames)
netlist_fnames<-sprintf("%s/%s/bloblist_DJF",input_directory,Varnames)
for (i in 1:length(stitchlist_fnames)){
  system(sprintf("ls %s/%s/%s > %s",input_directory,Varnames[i],stitch_searchstring,stitchlist_fnames[i]))
  system(sprintf("ls %s/%s/%s > %s",input_directory,Varnames[i],nostitch_searchstring,nostitch_fnames[i]))
  if (Varnames[i]!="MERRA" & Varnames[i]!="ERA" & Varnames[i]!="JRA" & Varnames[i]!="CFSR"){
    system(sprintf("ls %s/%s/%s > %s",input_directory,Varnames[i],netcdf_searchstring,netlist_fnames[i]))    
  }else{
    system(sprintf("ls %s/%s/%s > %s",input_directory,Varnames[i],netcdf2_searchstring,netlist_fnames[i]))
  }

}
##########################

#List of BlobStats files (from StitchBlobs output) to read into readfiles.R
stitch_lists<-stitchlist_fnames
#List of BlobStats files (from DetectBlobs output) to read into readfiles.R
detect_lists<-nostitch_fnames
#List of StitchBlobs files to read into readnetcdf.R
stitchblob_lists<-netlist_fnames
#Name of the StitchBlobs variable in the NetCDF file
varvec<-rep("Z_BLOB",length(Varnames))
#Name of the time, lat, lon axes
timename<-"time"
latname<-"lat"
lonname<-"lon"
#Transform the lon axis?
#from 0/360 to -180/180
transformto180<-FALSE
#from -180/180 to 0/360
transformto360<-TRUE
#Subset lat and lon if desired
minlat<- 25
maxlat<- 75
minlon<- 130
maxlon<- 270
#Regrid to 1 degree?
regridto1degree<-TRUE

