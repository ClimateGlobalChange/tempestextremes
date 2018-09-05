#NAMELIST FOR CREATING A SUMMARY TABLE USING --SUMMARIZE

#This is an example namelist to provide various arguments for summarizing
# blob information using --summarize
#Please do not change the variable names!!!

#########################
###FILE SPECIFICATIONS###
#########################
nrun_st<-1
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
#Options: Rdata ("R"), text ("text"), CSV ("CSV") or "session"
# Note that "session" can only be used in conjunction with --readfiles or --readtable
ftype_st="R"

#Leave a blank string if not reading in that particular file input type
#SINGLE FILE
filename_summ<-paste(input_dir,"table_merged.RData",sep="/")

#FILE LIST
filelist_summ<-""

#Keep or omit blobs that are comprised of multiple blobs?
#Default is TRUE
keepmerge=TRUE

#OUTPUT FILE NAME--leave blank if you don't want an output file
#RData file type-- will save a data.frame object
rfn_summ<-paste(output_dir,"table_summ.RData",sep="/")
#Name of the data frame variable (MANDATORY)
df_summ<-"df_summ_ERA_MERRA"

#text file type-- 
txt_summ<-""

#CSV
csv_summ<-""