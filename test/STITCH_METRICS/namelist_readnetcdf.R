#NAMELIST FOR READING IN NETCDF DATA USING --READNETCDF

#This is an example namelist to provide various arguments for reading NetCDF
# data into R using --readnetcdf
#Please do not change the variable names!!!
nrun_rn<-2
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
filename_netcdf<-c("~/tempestextremes/test/STITCH_METRICS/ERA_1980_DJF_comb_Z_blobs.nc",
                   "~/tempestextremes/test/STITCH_METRICS/MERRA_1980_DJF_comb_Z_blobs.nc")

#FILE LIST
filelist_netcdf<-""

#List of variable names to be read into R
#This is a vector in the form c("VAR1","VAR2","VAR3") etc
varvec<-"Z_BLOB"

#Optional list of output variable names
#If desired to leave it the same, just set to varvec
# outvec<-varvec
outvec<-list("ERA_BLOB","MERRA_BLOB")

#OUTPUT FILE NAME--leave blank if you don't want an output file
#RData file object-- will save an R workspace with 
# the NetCDF axes as vectors and the variables as 3D or 4D arrays
outrdata<-c("ERA_DJF.RData","MERRA_DJF.RData")
#NetCDF file
outnetcdf<-""

#Name of the NetCDF axes
timename<-"time"
levname<-"lev"
latname<-"lat"
lonname<-"lon"

#Convert the longitude axis range?
#Range from -180 to 180
transformto180<-FALSE
#Range from 0 to 360
transformto360<-TRUE

#Optional subsetting data
#If not subsetting, set the variable to ""
#Latitude axis
minlat<-25
maxlat<-75
#Longitude axis
minlon<-130
maxlon<-270
#Vertical level axis
minlev<-""
maxlev<-""