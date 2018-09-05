#NAMELIST FOR COMBINING DATA FROM STITCHBLOBS AND DETECTBLOBS USING --MERGETABLE

#This is an example namelist to provide various arguments for reading in data
# from StitchBlobs and DetectBlobs data
#Please do not change the variable names!!!
nrun_mt<-1
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
ftype_mt="R"

#Single file inputs
stitch_file<-""
detect_file<-""

#File lists
stitch_list<-""
detect_list<-""

#OUTPUT FILE NAME--leave blank if you don't want an output file
#RData file type-- will save a data.frame object
rfn_merged<-paste(output_dir,"table_merged.RData",sep="/")
#Name of the data frame variable (MANDATORY)
df_merged<-"df_merged_ERA_MERRA"

#text file type-- 
txt_merged<-""

#CSV
csv_merged<-""