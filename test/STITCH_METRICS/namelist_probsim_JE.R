#NAMELIST FOR COMPARING BLOCKING CO-OCCURRENCE USING --INTERCOMPARISON

#This is an example namelist to provide various arguments for calculating probability
# of co-occurrence, spatial similarity and pearson correlation using --intercomparison
#Please do not change the variable names!!!
nrun_ic<-2
#########################
###FILE SPECIFICATIONS###
#########################
#Location of working directory
work_dir<-"~/tempestextremes/test/STITCH_METRICS"
#Location where the input data files are stored
input_dir<-"~/BLOBSTATS_FILES"
#Location where the output data files are stored
output_dir<-"~/BLOBSTATS_FILES"

##################################
###INPUT/OUTPUT SPECIFICATIONS###
#################################
#Leave a blank string if not reading in that particular file input type
#NOTE: THE COARSER RESOLUTION FILE MUST BE DATASET 1
subdirs<-c("JRA","ERA")
SEASON<-c("DJF","JJA")
SECTOR<-c("NA","SP")
#Output from --mergetable or --readfiles containing the formatted per-timestep blocking data
table_file_1<-sprintf("%s/JRA/JRA_%s_%s_merged_table.RData",input_dir,SEASON,SECTOR)
#Name of the data frame
df_name_1<-"df_merged"
#If all of your data is contained in a single file, leave this as a blank string
table_file_2<-sprintf("%s/ERA/ERA_%s_%s_merged_table.RData",input_dir,SEASON,SECTOR)
df_name_2<-"df_merged"

#Blob data-- must be in RData format!
#Run the --readnetcdf function first
#If regrid is set to TRUE, regridding will be run
regrid=TRUE
#Name of the file for dataset 1 (MUST BE COARSER RESOLUTION, IF APPLICABLE)
blob_file_1<-sprintf("%s/JRA/JRA_%s_%s_blobdata.RData",input_dir,SEASON,SECTOR)
#Name of the array variable for dataset 1
var_name_1<-"Z_BLOB"
#Name of the file for dataset 2
blob_file_2<-sprintf("%s/ERA/ERA_%s_%s_blobdata.RData",input_dir,SEASON,SECTOR)
#name of the array variable for dataset 2
var_name_2<-"Z_BLOB"

rfn_ps=sprintf("%s/%s_%s_ic_JRA_ERA.RData",output_dir,SEASON,SECTOR)
txt_overlaps=sprintf("%s/%s_%s_ic_JRA_ERA_table.txt",output_dir,SEASON,SECTOR)
txt_ps=sprintf("%s/%s_%s_ic_JRA_ERA_probsim.txt",output_dir,SEASON,SECTOR)


