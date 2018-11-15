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
metadata_datasets="This is the results summary for DJF NA (December 1980-February 2005). 
There are 4 reanalysis datasets (ERA-Interim, MERRA2, JRA, and CFSR) and 12 historical model outputs. 
All of the reanalysis data is 6-hourly, and all of the model data is daily (12Z). \n\n
Intercomparison analysis for probability of co-occurrence, spatial similarity, Pearson pattern correlation, 
and RMSE is performed using only the common subset of timesteps, and all data has
been regridded to 1x1 degree (will assess whether changing the regrid resolution affects 
intercomparison results).
The individual dataset information for number of blobs, number of blocked days,
etc retains all time steps."
#This is the number of years in your dataset!
nyears<-2004-1980+1

##########################################
#USER_DEFINED
#Search string for stitch files
stitch_searchstring<-"*DJF_NA*stats.txt"
nostitch_searchstring<-"*DJF_NA*stats_nostitch.txt"
netcdf_searchstring<-"*DJF_NA*blobs.nc"
netcdf2_searchstring<-"*DJF*blobs.nc"
#########################################

###########
#FILE INFO
#Will always output RData files
#But there is an option to also output text files
# or CSV files with the table data
#Set to TRUE if you wish to have one or both of these
output_txt<-FALSE
output_csv<-FALSE

#Use DetectBlobs inputs?
use_detectblob<-TRUE

#Names of the datasets, which will be used to differentiate between them
#in the tables
#Varnames<-c("MERRA","JRA")
Varnames<-c("CMCC-CESM","CMCC-CM","CNRM-CM5","CanESM2","GFDL-CM3","GFDL-ESM2M",
            "MIROC-ESM","MIROC5","MPI-ESM-MR","MRI-CGCM3",
            "MRI-ESM1","NorESM1-M","CFSR","ERA","MERRA","JRA")
resolutions<-c("3.3341x3.75","TM159 (0.75)","T85 (1.39)","T42 (2.79)",
               "2x2.5","2x2.5","T42 (2.8)","T85 (1.39)","T63 (1.88)",
               "T106 (1.12)","T106 (1.12)","1.9x2.5","0.5x0.5","1x1",
               "0.5x0.625","1.25x1.25head")
#Working directory (where all of the R functions are)
work_dir<-"~/tempestextremes/test/STITCH_METRICS"
#Input directory
input_dir<-"~/BLOBSTATS_FILES"
#Output directory (where all files will go)
output_dir<-"~/BLOBSTATS_FILES"
#Name of output subdirectory
#This will be created in the output directory
output_subdir<-"DJF_NA"
#Name of the file prefix 
output_prefix<-"DJF_NA"

#Output name for the master namelist that will be generated
fname_namelist<-"DJF_NA_namelist_master.R"

#NOTE THAT ALL OF THE FILES CONTAINED IN THESE LISTS SHOULD HAVE FULL PATHNAMES
#Note that the input directory is appended to these list files!


#List of BlobStats files (from StitchBlobs output) to read into readfiles.R
#USER_DEFINED
#Make the lists on the fly using the search strings above!
##########################
stitchlist_fnames<-sprintf("%s/%s_stitch_list_NA",Varnames,Varnames)
nostitch_fnames<-sprintf("%s/%s_nostitch_list_NA",Varnames,Varnames)
netlist_fnames<-sprintf("%s/bloblist_DJF",Varnames)
for (i in 1:length(stitchlist_fnames)){
  system(sprintf("ls %s/%s/%s > %s/%s",input_dir,Varnames[i],stitch_searchstring,output_dir,stitchlist_fnames[i]))
  system(sprintf("ls %s/%s/%s > %s/%s",input_dir,Varnames[i],nostitch_searchstring,output_dir,nostitch_fnames[i]))
  if (Varnames[i]!="MERRA" & Varnames[i]!="ERA" & Varnames[i]!="JRA" & Varnames[i]!="CFSR"){
    system(sprintf("ls %s/%s/%s > %s/%s",input_dir,Varnames[i],netcdf_searchstring,output_dir,netlist_fnames[i]))    
  }else{
    system(sprintf("ls %s/%s/%s > %s/%s",input_dir,Varnames[i],netcdf2_searchstring,output_dir,netlist_fnames[i]))
  }

}
##########################

#stitch_lists<-c("JRA/JRA_stitch_list_NA","ERA/ERA_stitch_list_NA",
#                "MERRA/MERRA_stitch_list_NA","CFSR/CFSR_stitch_list_NA")
stitch_lists<-stitchlist_fnames
#List of BlobStats files (from DetectBlobs output) to read into readfiles.R
#If you are not using BlobStats files from DetectBlobs, set this variable to ""
#detect_lists<-c("JRA/JRA_nostitch_list_NA","ERA/ERA_nostitch_list_NA",
#                "MERRA/MERRA_nostitch_list_NA","CFSR/CFSR_nostitch_list_NA")
detect_lists<-nostitch_fnames
#List of StitchBlobs files to read into readnetcdf.R
#stitchblob_lists<-c("JRA/bloblist_DJF","ERA/bloblist_DJF",
#                    "MERRA/bloblist_DJF","CFSR/bloblist_DJF")
stitchblob_lists<-netlist_fnames

#Name of the StitchBlobs variable in the NetCDF file
varvec<-rep("Z_BLOB",length(Varnames))
#Name of the time, lat, lon axes
timename<-"time"
latname<-"lat"
lonname<-"lon"
#Transform the lon axis?
#from 0/360 to -180/180
transformto180<-TRUE
#from -180/180 to 0/360
transformto360<-FALSE
#Subset lat and lon if desired
minlat<-25
maxlat<-75
minlon<- -110
maxlon<-50
#Regrid to 1 degree?
regridto1degree<-TRUE

