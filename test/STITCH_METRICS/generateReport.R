require(knitr)
require(markdown)
require(rmarkdown)
require(reshape2)
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

#source("~/tempestextremes/test/STITCH_METRICS/namelist_report_JJA_SP.R")
#Generate the title string based on the variables

title_string<-paste("Comparison of blocking data for ",Varnames[1])
for (t in 2:length(Varnames)){
  title_string<-paste(title_string, Varnames[t], sep=", ")
}


md_file<-"report_template.Rmd"

avgdata<-data.frame(x=numeric(),y=numeric(),value=numeric(),
                     VAR=character(),lon=numeric(),lat=numeric())

#Make a list object that will have all of the relevant data
comparison_data<-list()
for (i in 1:length(Varnames)){
  comparison_data$varname[i]<-Varnames[i]
  #Load the merged table
  load(mergefiles[i])
  merge_dfname<-sprintf("V%d_merge",i)
  assign(merge_dfname,get(df_name))
  comparison_data$mergename[i]<-merge_dfname
  #load the summary table
  load(summfiles[i])
  summ_dfname<-sprintf("V%d_summ",i)
  assign(summ_dfname,df_summ)
  comparison_data$summname[i]<-summ_dfname
  #load the blob data
  load(blobfiles[i])
  assign(sprintf("lat%d",i),lat_axis)
  assign(sprintf("lon%d",i),lon_axis)
  assign(sprintf("time%d",i),time_format)
  assign(sprintf("blob%d",i),get(blobname[i]))
  temp_var<-get(sprintf("blob%d",i))
  temp_var[which(temp_var>0)]<-1
  assign(sprintf("blob%d",i),temp_var)
  #average the blob data
  avgname<-sprintf("avgblob%d",i)
  ablob<-apply(get(sprintf("blob%d",i)),c(1,2),mean)
  #Add to the long table for plotting
  temp<-melt(ablob,varnames=c("x","y"))
  temp$VAR<-rep(Varnames[i],nrow(temp))
  temp$lon<-lon_axis[temp$x]
  temp$lat<-lat_axis[temp$y]
  avgdata<-rbind(avgdata,temp)
  assign(avgname,ablob)
}

avgdata$VAR<-factor(avgdata$VAR,levels=Varnames)

#Generate the report from the template
rmarkdown::render(md_file,output_file=output_name)

