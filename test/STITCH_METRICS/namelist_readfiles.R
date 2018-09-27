#NAMELIST FOR READING IN DATA USING --READFILES

#This is an example namelist to provide various arguments for reading in data
#Please do not change the variable names!!!

nrun_rf<-1
#########################
###FILE SPECIFICATIONS###
#########################
#Location of working directory
work_dir<-"~/tempestextremes/test/STITCH_METRICS"
#Location where the input data files are stored
input_dir<-"~/tempestextremes/test/STITCH_METRICS"
#Location where the output data files are stored
output_dir<-"~/tempestextremes/test/STITCH_METRICS"

multiple_readfiles<-FALSE
##################################
###INPUT/OUTPUT SPECIFICATIONS###
#################################

#The string that will go in the "var" column
#Use this to distinguish between different datasets
varname<-"ERA"

#Leave a blank string if not reading in that particular file input type
#SINGLE FILE
#Input BlobStats file 
filename_stitchblobs<-""

#FILE LIST
#Input list of BlobStats files 
filelist_stitchblobs<-paste(input_dir,"stitch_list",sep="/")

#OUTPUT FILE NAME--leave blank if you don't want an output file
#RData file type-- will save a data.frame object
rfn_stitch<-paste(output_dir,"table_stitch.RData",sep="/")
#Name of the data frame variable
df_stitchname<-"df_stitch"

#text file type-- 
txt_stitch<-""

#CSV
csv_stitch<-""

