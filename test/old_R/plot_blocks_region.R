library(maps)
library(maptools)
for (region in c("
                 warNC","NP","SA","SI","SP")){
#for (region in c("NA","NC","NP","SA","SI","SP")){
  load(sprintf('~/block_r_data/%s_pv_z_ghg_block_data.RData',region))
  time_format1<-time_format
  load(sprintf('~/block_r_data/%s_pv_z_inst_data.RData',region))
  #Need to reverse the latitude axis for each of the data
  lats_seq<-lats_seq[length(lats_seq):1]
  z_inst<-z_inst[,length(lats_seq):1,]
  z_anom<-z_anom[,length(lats_seq):1,]
  pv_anom<-pv_anom[,length(lats_seq):1,]
  ghg<-ghg[,length(lats_seq):1,]
  
  if (region=="NA" | region=="SA"){
    lons_plot=lons_seq_c
    mapr="world"
  }else{
    lons_plot=lons_seq
    mapr="world2"
  }
  
  #NOTE!!!
  #the time vectors for blobs and z/pv are different!!
  #blobs starts in March 1980
  #z/pv starts in January 1980
  
  months<-format(time_format,"%b")
  years<-format(time_format,"%Y")
  
  
  months2<-format(time_format2,"%b")
  years2<-format(time_format2,"%Y")
  
  for (season in c("MAM","JJA","SON","DJF")){
    load(sprintf('~/block_r_data/stats_stitch_%s_%s_table.RData',season,region))
    
    dir.create(file.path("~/pics_test_case/", sprintf("%s_%s",region,season)), showWarnings = FALSE)

    #Pick out some specific cases of blocking that occurred! Use the stats table 
    # for reference dates
    
    uniq_dates<-unique(df_tot_stitch$date)
    
    hgt.cols<-colorRampPalette(c("purple","blue","cyan4","green","yellow", "orange","red","darkred"))(34)
    for (u in 1:length(uniq_dates)){
      uniq=uniq_dates[u]
      check_date<-format(uniq,"%m_%d")
      if (check_date != "02_29"){
        t1=which(time_format==uniq)
        t2=which(time_format2==uniq)
        if (length(t1)>0 & length(t2)>0){
          fname<-sprintf("~/pics_test_case/%s_%s/%s_%s_%02dZ_noT.png",region,season,region,time_format[t1],time_hrs[t1])
          print(fname)
          png(fname,height=600,width=800)
          
          map(mapr,xlim=range(lons_plot),ylim=range(lats_seq),fill=TRUE,col="grey")
          title(sprintf("Z500 %s %02dZ, PV* (green) Z* (blue) ZG (purple)",time_format[t1],time_hrs[t1]))
          map.axes()
          contour(lons_plot,lats_seq,z_inst[,,t2],levels=seq(4500,6100,50),drawlabels=TRUE,add=TRUE,col=hgt.cols,lwd=2)
          contour(lons_plot,lats_seq,pv_anom[,,t1],levels=c(0,1),add=TRUE,col="chartreuse4",drawlabels=FALSE,lwd=5)
          contour(lons_plot,lats_seq,z_anom[,,t1],levels=c(0,1),add=TRUE,col="cornflowerblue",drawlabels=FALSE,lwd=5)
          contour(lons_plot,lats_seq,ghg[,,t1],levels=c(0,1),add=TRUE,col="purple",drawlabels=FALSE,lwd=5)
          dev.off()
        }

      }

    }

  }
}