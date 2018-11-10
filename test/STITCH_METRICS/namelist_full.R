#This is an example namelist for an output report
#The requirements are as follows:

#1) StitchBlobs files from each dataset
#2) Corresponding BlobStats files from each dataset
#3) Optional but recommended: BlobStats output from
# the separate DetectBlobs function
# This will help deal with cases where two blobs
# merge into a single blob at a later time 

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
Varnames<-c("JRA","ERA","MERRA","CFSR")

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
stitch_lists<-c("JRA/JRA_stitch_list_NA","ERA/ERA_stitch_list_NA",
                "MERRA/MERRA_stitch_list_NA","CFSR/CFSR_stitch_list_NA")

#List of BlobStats files (from DetectBlobs output) to read into readfiles.R
#If you are not using BlobStats files from DetectBlobs, set this variable to ""
detect_lists<-c("JRA/JRA_nostitch_list_NA","ERA/ERA_nostitch_list_NA",
                "MERRA/MERRA_nostitch_list_NA","CFSR/CFSR_nostitch_list_NA")

#List of StitchBlobs files to read into readnetcdf.R
stitchblob_lists<-c("JRA/bloblist_DJF","ERA/bloblist_DJF",
                    "MERRA/bloblist_DJF","CFSR/bloblist_DJF")
#Name of the StitchBlobs variable in the NetCDF file
varvec<-rep("Z_BLOB",4)
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


