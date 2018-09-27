#THIS IS THE MASTER NAMELIST FROM ALL OF THE VARIOUS FUNCTIONS
#SEE EACH INDIVIDUAL NAMELIST FOR DESCRIPTIONS

#########################
###FILE SPECIFICATIONS###
#########################
#Location of working directory
work_dir<-"~/tempestextremes/test/STITCH_METRICS"
#Location where the input data files are stored
input_dir<-"/group/paullricgrp2/ERA_MCP/ERA_BLOB"
#Location where the output data files are stored
output_dir<-"/group/paullricgrp2/ERA_MCP/ERA_BLOB"

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE THE READFILES OPERATION IS RUN 4 TIMES
#The analysis is broken up by region because otherwise it will take too long!
SECTOR=c("NA", "NC", "NP", "SA", "SI" ,"SP")
#READFILES-------------
nrun_rf<-12
###USER-DEFINED###
#This will be run 12 times 
stitch_search_str<-sprintf("%s/ERA*%s_Z_stats.txt",input_dir,SECTOR)
nostitch_search_str<-sprintf("%s/ERA*%s_Z_stats_nostitch.txt",input_dir,SECTOR)
search_str<-c(stitch_search_str,nostitch_search_str)

list_files<-c(sprintf("%s/ERA_stitch_list_%s",input_dir,SECTOR),sprintf("%s/ERA_nostitch_list_%s",input_dir,SECTOR))
list_rnames<-c(sprintf("%s/ERA_%s_stitchtable.Rdata",output_dir,SECTOR),sprintf("%s/ERA_%s_nostitchtable.Rdata",output_dir,SECTOR))
list_dfnames<-"table_out"
for (i in length(list_files)){
  system(sprintf("ls %s > %s",search_str[i],list_files[i]))
}
#################

varname<-rep("ERA",12)
filename_stitchblobs<-""
filelist_stitchblobs<-list_files
rfn_stitch<-list_rnames
df_stitchname<-list_dfnames
txt_stitch<-""
csv_stitch<-""

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE THE READTABLE OPERATION COMBINES 
#THE STITCH AND THE NOSTITCH DATA INTO TWO DISTINCT TABLES
# THEREFORE THE OPERATION IS RUN TWICE

#READTABLE-------------
nrun_rt<-2
#USER-DEFINED
#Write the file names above to lists for the StitchBlobs and DetectBlobs outputs
#File names
stitch_rlist<-paste(input_dir,"stitch_files_list",sep="/")
nostitch_rlist<-paste(input_dir,"nostitch_files_list",sep="/")
#Write to file
writeLines(rfn_stitch[c(1,3)],stitch_rlist)
writeLines(rfn_stitch[c(2,4)],nostitch_rlist)
#Names of output files
list_out<-c("table_stitch.RData","table_nostitch.RData")
############

ftype_rt="R"
filename_read<-""
filelist_read<-c(stitch_rlist,nostitch_rlist)
rfn_combine<-paste(output_dir,list_out,sep="/")
df_combinename<-c("df_stitch_ERA_MERRA","df_nostitch_ERA_MERRA")
txt_combine<-""
csv_combine<-""

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE MERGETABLE IS USED
#TO COMBINE STITCHBLOBS AND DETECTBLOBS OUTPUTS FROM BOTH ERA AND MERRA
#INTO ONE FILE
#NOTE THAT THE FILE NAMES FROM READTABLE WERE USED 
#TO DEFINE THE VARIABLES HERE

#MERGETABLE------------
nrun_mt<-1
ftype_mt="R"
stitch_file<-rfn_combine[1]
detect_file<-rfn_combine[2]
stitch_list<-""
detect_list<-""
rfn_merged<-paste(output_dir,"table_merged.RData",sep="/")
df_merged<-"df_merged_ERA_MERRA"
txt_merged<-""
csv_merged<-""

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE SUMMARIZE IS USED
# TO SUMMARIZE THE MERGETABLE OUTPUT FROM THE PREVIOUS SECTION
#SUMMARIZE----------
nrun_st<-1
ftype_st<-"R"
filename_summ<-rfn_merged
filelist_summ<-""
keepmerge<-TRUE
rfn_summ<-paste(output_dir,"table_summ.RData",sep="/")
df_summ<-"df_summ_ERA_MERRA"
txt_summ<-""
csv_summ<-""

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE NETCDF DATA
#IS READ INTO AN ARRAY 
#The min and max boundaries are defined per region below

#READNETCDF-------------
nrun_rn<-3
########################
#USER DEFINED#



###################################################################

filename_netcdf<-""
filelist_netcdf<-""
varvec<-list(blob_vars,blob_vars,wind_vars)
outvec<-varvec
outrdata<-c("blob_data_MERRA.RData","blob_data_ERA.RData","wind_data.RData")
outnetcdf<-""
timename<-"time"
levname<-"lev"
latname<-"lat"
lonname<-"lon"
#Since this is in the Pacific, we want the longitude range to be in the 0 to 360 range
#MERRA has a longitude axis from -180 to 180 while ERA has a longitude axis from 0 to 360
transformto180<-FALSE
transformto360<-TRUE
minlat<-25
maxlat<-75
minlon<-130
maxlon<-270
minlev<-c("","",500)
maxlev<-c("","",500)