#THIS IS THE MASTER NAMELIST FROM ALL OF THE VARIOUS FUNCTIONS
#SEE EACH INDIVIDUAL NAMELIST FOR DESCRIPTIONS

#########################
###FILE SPECIFICATIONS###
#########################
#Location of working directory
work_dir<-"~/tempestextremes/test/STITCH_METRICS"
#Location where the input data files are stored
input_dir<-"~/tempestextremes/test/STITCH_METRICS"
#Location where the output data files are stored
output_dir<-"~/tempestextremes/test/STITCH_METRICS"

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE THE READFILES OPERATION IS RUN 4 TIMES

#READFILES-------------
nrun_rf<-4
###USER-DEFINED###
#This will be run 4 times 
list_files<-c("stitch_list","nostitch_list","merra_list","merra_nolist")
list_rnames<-c("table_estitch.RData","table_enostitch.RData","table_mstitch.RData","table_mnostitch.RData")
list_dfnames<-c("df_estitch","df_enostitch","df_mstitch","df_mnostitch")
#################

nhrs<-6
varname<-c("ERA","ERA","MERRA","MERRA")
filename_stitchblobs<-""
filelist_stitchblobs<-paste(input_dir,list_files,sep="/")
rfn_stitch<-paste(output_dir,list_rnames,sep="/")
df_stitchname<-list_dfnames
txt_stitch<-""
csv_stitch<-""

#THIS IS AN EXAMPLE OF A NAMELIST SECTION WHERE THE READTABLE OPERATION COMBINES 
#THE STITCH AND THE NOSTITCH DATA INTO TWO DISTINCT TABLES
# THEREFORE THE OPERATION IS RUN TWICE

#READTABLE-------------
nrun_rt<-2
#USER-DEFINED
list_files<-c("table_stitch_in","table_nostitch_in")
list_out<-c("table_stitch.RData","table_nostitch.RData")
############

ftype_rt="R"
filename_read<-""
filelist_read<-paste(input_dir,list_files,sep="/")
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
#IS READ INTO AN ARRAY AND SUBSET TO ONLY THE NORTH PACIFIC
#THE FIRST TWO RUNS READ IN STITCHBLOBS DATA FROM ERA AND MERRA
#THE SECOND RUN READS IN A LIST OF FILES CONTAINING WIND DATA 
#FROM ERA FOR JJA
#NOTE THAT WE USE A LIST TO COMBINE VECTORS OF VARIABLE NAMES!!

#READNETCDF-------------
nrun_rn<-3
########################
#USER DEFINED#
flist_windfiles<-paste(input_dir,"ERA_JJA_list",sep="/")
blobfile_m<-paste(input_dir,"MERRA_1980_JJA_comb_Z_blobs.nc",sep="/")
blobfile_e<-paste(input_dir,"ERA_1980_JJA_comb_Z_blobs.nc",sep="/")
blob_vars<-c("Z_BLOB")
wind_vars<-c("U","V")
###################################################################

filename_netcdf<-c(blobfile_m,blobfile_e,"")
filelist_netcdf<-c("","",flist_windfiles)
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