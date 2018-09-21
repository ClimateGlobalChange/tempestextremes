require(knitr)
require(markdown)
require(rmarkdown)

#Returns in the -180->180 range
lon_convert<-function(lon){
  distFrom180=lon-180.
  return(ifelse(
    distFrom180<0,
    lon,
    -(180-distFrom180)
  ))
}
#Returns in the 0->360 range
lon_convert2<-function(lon){
  return(ifelse(lon<0,360+lon,lon))
}

#Load the two datasets
V1<-"ERA"
V2<-"MERRA"

title_string<-sprintf("Comparison of blocking data from %s and %s",V1,V2)
md_file<-"~/tempestextremes/test/STITCH_METRICS/report_template.Rmd"

#Data tables for V1 and V2
load("table_merged.RData")
#This is the per-timestep data
df_table<-get(df_name)
#This is the summarized data
#MAKE A NOTE TO ADD THE DF NAME TO THE R FILE
load("table_summ.RData")
df_summ<-df_summ_ERA_MERRA

#Load the blob data
load("ERA_DJF.RData")
ERA_lat<-lat_axis
ERA_lon<-lon_axis
ERA_BLOB[which(ERA_BLOB>0)]<-1
V1_dens<-apply(ERA_BLOB,c(1,2),mean)
load("MERRA_DJF.RData")
MERRA_lat<-lat_axis
MERRA_lon<-lon_axis
MERRA_BLOB[which(MERRA_BLOB>0)]<-1
V2_dens<-apply(MERRA_BLOB,c(1,2),mean)

#Generate the report from the template
rmarkdown::render(md_file,output_file="test.html")

