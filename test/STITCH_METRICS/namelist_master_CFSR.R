#THIS IS THE MASTER NAMELIST FROM ALL OF THE VARIOUS FUNCTIONS
#SEE EACH INDIVIDUAL NAMELIST FOR DESCRIPTIONS

#########################
###FILE SPECIFICATIONS###
#########################
#Location of working directory
work_dir<-"~/tempestextremes/test/STITCH_METRICS"
#Location where the input data files are stored
input_dir<-"~/BLOBSTATS_FILES/CFSR"
#Location where the output data files are stored
output_dir<-"~/BLOBSTATS_FILES/CFSR"

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE THE READFILES OPERATION IS RUN 4 TIMES
#The analysis is broken up by region because otherwise it will take too long!
#We're looking at Northern Hemisphere Pacific and Southern Hemisphere Pacific
#Winter is DJF in NH and JJA in SH
SECTOR=c("NA", "SP")
SEASON=c("DJF","JJA")
#READFILES-------------
nrun_rf<-4
###USER-DEFINED###
#This will be run 4 times 
stitch_search_str<-sprintf("%s/CFSR*_%s_%s_Z_stats.txt",input_dir,SEASON,SECTOR)
nostitch_search_str<-sprintf("%s/CFSR*%s_%s_Z_stats_nostitch.txt",input_dir,SEASON,SECTOR)
search_str<-c(stitch_search_str,nostitch_search_str)

list_files<-c(sprintf("%s/CFSR_stitch_list_%s",input_dir,SECTOR),sprintf("%s/CFSR_nostitch_list_%s",input_dir,SECTOR))
list_rnames<-c(sprintf("%s/CFSR_%s_%s_stitchtable.Rdata",output_dir,SEASON,SECTOR),
               sprintf("%s/CFSR_%s_%s_nostitchtable.Rdata",output_dir,SEASON,SECTOR))
list_dfnames<-"table_out"
for (i in 1:length(list_files)){
  system(sprintf("ls %s > %s",search_str[i],list_files[i]))
}
#################

varname<-rep("CFSR",4)
filename_stitchblobs<-""
filelist_stitchblobs<-list_files
rfn_stitch<-list_rnames
df_stitchname<-list_dfnames
txt_stitch<-""
csv_stitch<-""

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE THE READTABLE OPERATION COMBINES 
#THE STITCH AND THE NOSTITCH DATA INTO TWO DISTINCT TABLES
# THEREFORE THE OPERATION IS RUN TWICE
#(THIS IS OLD INPUT AND DOES NOT MATCH UP WITH PREVIOUS SECTION)

#READTABLE-------------
nrun_rt<-2
#USER-DEFINED
#Write the file names above to lists for the StitchBlobs and DetectBlobs outputs
#File names
# stitch_rlist<-paste(input_dir,"stitch_files_list",sep="/")
# nostitch_rlist<-paste(input_dir,"nostitch_files_list",sep="/")
# #Write to file
# writeLines(rfn_stitch[c(1,3)],stitch_rlist)
# writeLines(rfn_stitch[c(2,4)],nostitch_rlist)
# #Names of output files
# list_out<-c("table_stitch.RData","table_nostitch.RData")
# ############
# 
# ftype_rt="R"
# filename_read<-""
# filelist_read<-c(stitch_rlist,nostitch_rlist)
# rfn_combine<-paste(output_dir,list_out,sep="/")
# df_combinename<-c("df_stitch_ERA_MERRA","df_nostitch_ERA_MERRA")
# txt_combine<-""
# csv_combine<-""

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE MERGETABLE IS USED
#TO COMBINE STITCHBLOBS AND DETECTBLOBS OUTPUTS 

#USER-DEFINED
#WE WILL BE WRITING TWO MERGED FILES FOR EACH OF THE TWO BASINS


#MERGETABLE------------
nrun_mt<-2
ftype_mt="R"
stitch_file<-list_rnames[1:2]
detect_file<-list_rnames[3:4]
stitch_list<-""
detect_list<-""
rfn_merged<-sprintf("%s/CFSR_%s_%s_merged_table.RData",output_dir,SEASON,SECTOR)
df_merged_name<-"df_merged"
txt_merged<-""
csv_merged<-""

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE SUMMARIZE IS USED
# TO SUMMARIZE THE MERGETABLE OUTPUT FROM THE PREVIOUS SECTION
#SUMMARIZE----------
nrun_st<-2
ftype_st<-"R"
filename_summ<-rfn_merged
filelist_summ<-""
keepmerge<-TRUE
rfn_summ<-sprintf("%s/CFSR_%s_%s_summ_table.RData",output_dir,SEASON,SECTOR)
df_summ_name<-"df_summ"
txt_summ<-""
csv_summ<-""

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE NETCDF DATA
#IS READ INTO AN ARRAY 
#The min and max boundaries are defined per region below

#READNETCDF-------------
nrun_rn<-2
########################
#USER DEFINED#
netcdf_search_str<-sprintf("%s/CFSR*_%s_comb_Z_blobs.nc",input_dir,SEASON)
list_netcdf_files<-sprintf("%s/bloblist_%s",input_dir,SEASON)
for (i in 1:length(list_netcdf_files)){
  system(sprintf("ls %s > %s",netcdf_search_str[i],list_netcdf_files[i]))
}
rfiles_blobdata<-sprintf("%s/CFSR_%s_%s_blobdata.RData",output_dir,SEASON,SECTOR)
########################

filename_netcdf<-""
filelist_netcdf<-list_netcdf_files
varvec<-c("Z_BLOB")
outvec<-varvec
outrdata<-rfiles_blobdata
outnetcdf<-""
timename<-"time"
levname<-"lev"
latname<-"lat"
lonname<-"lon"
#Subsetting the regions/seasons in both instances
#For MERRA, need to change axis in order to match ERA
transformto180<-c(TRUE,FALSE)
transformto360<-c(FALSE,TRUE)
minlat<-c(25,-75)
maxlat<-c(75,-25)
minlon<-c(-110,120)
maxlon<-c(50,310)
minlev<-""
maxlev<-""

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE TWO DATASETS ARE COMPARED
