#NAMELIST FOR READING IN EXISTING DATA USING --READTABLE

#This is an example namelist to provide various arguments for reading in data
# from existing --readfiles output
#Please do not change the variable names!!!
nrun_rt<-1
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

#The file type that is being read in
#Options: Rdata ("R"), text ("text") or CSV ("CSV")
ftype_rt="R"

#Leave a blank string if not reading in that particular file input type
#SINGLE FILE
#Input BlobStats file from StitchBlobs (--readfiles)
filename_read<-paste(input_dir,"table_stitch.RData",sep="/")

#FILE LIST
#Input list of BlobStats files from StitchBlobs (--readfiles)
filelist_read<-""

#OUTPUT FILE NAME--leave blank if you don't want an output file
#RData file type-- will save a data.frame object
rfn_combine<-paste(output_dir,"table_combine.RData",sep="/")
#Name of the data frame variable (MANDATORY)
df_combinename<-"df_stitch_ERA_MERRA"

#text file type-- 
txt_combine<-""

#CSV
csv_combine<-""