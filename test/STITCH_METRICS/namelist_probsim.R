#NAMELIST FOR COMPARING BLOCKING CO-OCCURRENCE USING --PROBSIM

#This is an example namelist to provide various arguments for calculating probability
# of co-occurrence and spatial similarity using --probsim
#Please do not change the variable names!!!
nrun_ps<-1
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
#NOTE: THE COARSER RESOLUTION FILE MUST BE DATASET 1

#Output from --readfiles containing the formatted per-timestep blocking data
table_file_1<-paste(input_dir,"table_summ.RData",sep="/")
#Name of the data frame
df_name_1<-"df_summ_ERA_MERRA"
#If all of your data is contained in a single file, leave this as a blank string
table_file_2<-""
df_name_2<-""

#Blob data-- must be in RData format!
#Run the --readnetcdf function first
#Name of the file
blob_file_1<-paste(input_dir,"ERA_1980_DJF_comb_Z_blobs.nc",sep="/")
#Name of the array variable
var_name_1<-"ERA_BLOB"

blob_file_2<-paste(input_dir,"MERRA_1980_DJF_comb_Z_blobs.nc",sep="/")
var_name_2<-"MERRA_BLOB"



