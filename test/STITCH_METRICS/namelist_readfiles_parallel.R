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