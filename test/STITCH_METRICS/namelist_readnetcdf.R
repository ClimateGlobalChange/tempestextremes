#NAMELIST FOR READING IN NETCDF DATA USING --READNETCDF

#This is an example namelist to provide various arguments for reading NetCDF
# data into R using --readnetcdf
#Please do not change the variable names!!!
nrun_rn<-1
#########################
###FILE SPECIFICATIONS###
#########################
#Location of working directory
work_dir<-"~/tempestextremes/test/STITCH_METRICS"
#Location where the input data files are stored
input_dir<-"~/tempestextremes/test/STITCH_METRICS"
#Location where the output data files are stored
output_dir<-"~/tempestextremes/test/STITCH_METRICS"

##################################
###INPUT/OUTPUT SPECIFICATIONS###
#################################
#Leave a blank string if not reading in that particular file input type
#SINGLE FILE
filename_netcdf<-paste(input_dir,"ERA_1980_JJA_comb_Z_blobs.nc",sep="/")

#FILE LIST
filelist_netcdf<-""

#List of variable names to be read into R
#This is a vector in the form c("VAR1","VAR2","VAR3") etc
varvec<-c("Z_BLOB")

#Optional list of output variable names
#If desired to leave it the same, just set to varvec
# outvec<-varvec
outvec<-varvec

#OUTPUT FILE NAME--leave blank if you don't want an output file
#RData file object-- will save an R workspace with 
# the NetCDF axes as vectors and the variables as 3D or 4D arrays
outrdata<-"blob_data.RData"
#NetCDF file
outnetcdf<-""

#Name of the NetCDF axes
timename<-"time"
levname<-"lev"
latname<-"lat"
lonname<-"lon"

#Optional subsetting data
#If the entire latitude/longitude extent is desired, set these variables to NULL
#Example:
# minlat<-NULL

#Latitude axis
minlat<-25
maxlat<-75
#Longitude axis
minlon<-130
maxlon<-270
#Vertical level axis
minlev<-NULL
maxlev<-NULL