#The requirements are as follows:
#1) StitchBlobs files from each dataset
#2) Corresponding BlobStats files from each dataset
#3) Optional but recommended: BlobStats output from 
# the separate DetectBlobs function 
# (This will help deal with cases where two blobs 
# merge into a single blob at a later time)  

#REPORT HEADER: ANY TEXT THAT YOU WANT TO GO AT THE TOP OF THE REPORT 
#NOTE: must be in quotes. Careful of any special characters 
metadata_datasets="This is the results summary for DJF NP. Four of the datasets are reanalyses (ERA-Interim, CFSR, MERRA-2, JRA). Two are models: one is a historical dataset and one follows the RCP8.5 scenario.

All of the reanalysis is 6-hourly and all of the model data is daily (12Z). All data has been regridded to 1 degree in order to be able to do the spatial comparisons." 
#This is the number of years in each dataset! (Vector of integers) 
nyears<-c((2099-2006+1),25, 25, 25, 25, 25)

########### 
#FILE INFO 
#Will always output RData files 
#But there is an option to also output text files 
# or CSV files with the table data 
#Set following variables to TRUE if you wish to have one or both of these outputs 
output_txt<-FALSE 
output_csv<-FALSE  

#Names of the datasets (vector of strings) 
#Example: Varnames<-c("ERA-Interim","CFSR","MERRA-2") 
#For all of the vectors of strings, make sure that the lengths are identical to 
# the length of this vector! 
Varnames<-c("CMCC-CESM-RCP8.5","CMCC-CESM","JRA", "MERRA","CFSR","ERA")  

#The spatial resolutions of the datasets (Vector of strings.  
#Optional-- delete this variable if you don't want to use it) 
#Example: resolutions<-c("1x1","0.5x0.625","0.5x0.5") 
#resolutions<-c()  

#Working directory (String) 
#This is where all of the R function files are stored  
work_dir<-"~/tempestextremes/src/blobmetrics" 
#Output directory(String) 
# Main directory where all output files will go
output_dir<-"~/BLOBSTATS_FILES" 
#Name of output subdirectory (String) 
# This will be created in the output directory specified above 
output_subdir<-"MODEL_REAN" 
#Name of the file prefix (String).  
#Output file will be [prefix]_[Varname]_[suffix] (for example, [prefix]_ERA_stitchtable.RData)  
output_prefix<-"DJF_NP" 
#Output name for the master namelist that will be generated using this file template (String) 
fname_namelist<-"DJF_NP_model_reanalysis_dataset_namelist_master.R"  


#Input lists of files 
#Must use full pathname for each of these lists! 
#Example: stitch_lists<-c("~/input_dir/ERA/ERA_list","~/input_dir/CFSR/CFSR_list","~/input_dir/MERRA/MERRA_list") 
####USER DEFINED
stitch_searchstring<-"*DJF_NP*stats.txt"
nostitch_searchstring<-"*DJF_NP*stats_nostitch.txt"
netcdf_searchstring_model<-"*DJF_NP*blobs.nc"
netcdf_searchstring_rean<-"*DJF*blobs.nc"

input_directory<-"~/BLOBSTATS_FILES"
stitchlist_fnames<-sprintf("%s/%s/%s_stitch_list_DJF_NP",input_directory,Varnames,Varnames)
nostitchlist_fnames<-sub("stitch","nostitch",stitchlist_fnames)
netcdf_fnames<-sub("stitch","blob",stitchlist_fnames)
for (i in 1:length(stitchlist_fnames)){
  system(sprintf("ls %s/%s/%s>%s",input_directory,Varnames[i],
                 stitch_searchstring,stitchlist_fnames[i]))
  system(sprintf("ls %s/%s/%s>%s",input_directory, Varnames[i],
                 nostitch_searchstring,nostitchlist_fnames[i]))
  #I need to do this because I used a different filename convention for model vs reanalysis...
  if (i<3){
    system(sprintf("ls %s/%s/%s > %s",input_directory,Varnames[i],
                   netcdf_searchstring_model,netcdf_fnames[i]))
  }
  else{
    system(sprintf("ls %s/%s/%s > %s",input_directory, Varnames[i],
                   netcdf_searchstring_rean,netcdf_fnames[i]))
  }

}

##################
#List of BlobStats files (from StitchBlobs output) to read into readfiles.R (vector of strings) 
stitch_lists<-stitchlist_fnames
#Use DetectBlobs inputs? (TRUE/FALSE) 
use_detectblob<-TRUE  
#List of BlobStats files (from DetectBlobs output) to read into readfiles.R (vector of strings) 
#If not using BlobStats files from DetectBlobs, then delete this variable 
detect_lists<-nostitchlist_fnames
#List of StitchBlobs files to read into readnetcdf.R (vector of strings) 
stitchblob_lists<-netcdf_fnames
#Name of the StitchBlobs variable in the NetCDF file (vector of strings) 
#Note: if the variable name is identical for all file, you can do 
# varvec<-rep("VARNAME",length(Varnames)) 
varvec<-rep("Z_BLOB",length(Varnames))
#Name of the time, lat, lon axes (String) 
timename<-"time" 
latname<-"lat" 
lonname<-"lon" 
#Transform the lon axis? (TRUE/FALSE) 
#Note: if a dataset already has this longitude extent, it will do nothing 
#from 0/360 to -180/180 
transformto180<-FALSE 
#from -180/180 to 0/360 
transformto360<-TRUE 
#Subset lat and lon if desired (note: minlon and maxlon should correspond to 
# the appropriate longitude extent, i.e. either in the -180/180 range or 0/360 range) 
#If you don't wish to subset, delete these four variables 
minlat<- 25
maxlat<- 75 
minlon<- 130
maxlon<- 270 
#Regrid to 1 degree? (TRUE/FALSE) 
regridto1degree<-TRUE
#For Pearson correlation and RMSE: only use common subset of time?
# (Set to false if comparing two different time periods, such as historical vs RCP8.5)
useCommonTime<-FALSE

#Which sections will be included in the output report? (T/F)
#Initial summary table-- requires output from --summarize and --readfiles/--mergetable
includeSummTable<-TRUE
#Blocking frequency plots
includeFrequencyPlots<-TRUE
#Pearson Pattern Correlation between blocking frequencies
includePearson<-TRUE
#Root mean square error between blocking frequencies
includeRMSE<-TRUE
#Density plots for duration, speed, and size and associated p-values 
includeDensityP<-TRUE
#Intercomparison
includeProbability<-FALSE
includeSpatialSimilarity<-FALSE
