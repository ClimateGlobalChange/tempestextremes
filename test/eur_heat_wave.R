# source("~/tempestextremes/test/sector_funcs.R")
# 
# #First: anomalies per season per sector
# files_dir<-"/Volumes/ExFAT_drive/ERA_files"
# img_dir<-"/Users/mariellep/figs"
# 
# sectors<-c("NA","NC","NP","SA","SI","SP")
# colors<-c("red","orange","pink","blue","green","purple")
# 
# #datas<-c("climo","2xCO2","SSTplus2")
# LEFT_BOUND=c(250, 30, 130, 290, 20, 120)
# RIGHT_BOUND=c(50, 150, 270, 40, 140, 310)
# MIN_LAT=c(25, 25, 25, -75, -75, -75)
# MAX_LAT=c(75, 75, 75, -25, -25, -25)
# 
# f<-open.nc("/Volumes/ExFAT_drive/ERA_files/ERA_2003/ERA_2003_JJA_sfc.nc")
# lats<-read.nc(f)$latitude
# lons<-read.nc(f)$longitude
# timevar<-read.nc(f)$time
# 
# 
# #Get array indices for bounds
# li<-get_index(lons,LEFT_BOUND[1])
# ri<-get_index(lons,RIGHT_BOUND[1])
# bi<-get_index(lats,MIN_LAT[1])
# ti<-get_index(lats,MAX_LAT[1])
# lon_length<-ri+(length(lons)-li+1)
# lat_length<-(bi-ti)+1
# time_length<-length(timevar)
# 
# lon_seq<-seq(-110,50)
# lat_seq<-seq(25,75)
# 
# #Initialize array that will hold data
# t_sub<-array(NA,c(lon_length,lat_length,time_length))
# #Need to deal with fact that lat axis has decreasing values :(
# t_sub[1:(lon_length-ri),,]<-read.nc(f)$t2m[li:length(lons),bi:ti,]
# t_sub[(lon_length-ri+1):lon_length,,]<-read.nc(f)$t2m[1:ri,bi:ti,]
# 
# t_sub_c<-t_sub-273
# close.nc(f)
# 
# 
# pv_blob_sub<-array(NA,c(lon_length,lat_length,time_length))
# f<-open.nc(paste(files_dir,"/ERA_blobs/ERA_2003_JJA_NA_blobs.nc",sep=""))
# pv_blob_sub[1:(lon_length-ri),,]<-read.nc(f)$PV_BLOB[li:length(lons),bi:ti,]
# pv_blob_sub[(lon_length-ri+1):lon_length,,]<-read.nc(f)$PV_BLOB[1:ri,bi:ti,]
# close.nc(f)
# 
# 
# z_blob_sub<-array(NA,c(lon_length,lat_length,time_length))
# f<-open.nc(paste(files_dir,"/ERA_blobs/ERA_2003_JJA_NA_Zblobs.nc",sep=""))
# z_blob_sub[1:(lon_length-ri),,]<-read.nc(f)$Z_BLOB[li:length(lons),bi:ti,]
# z_blob_sub[(lon_length-ri+1):lon_length,,]<-read.nc(f)$Z_BLOB[1:ri,bi:ti,]
# close.nc(f)
# 
# gh_blob_sub<-array(NA,c(lon_length,lat_length,time_length))
# f<-open.nc(paste(files_dir,"/ERA_blobs/ERA_2003_JJA_NA_GHGblobs.nc",sep=""))
# gh_blob_sub[1:(lon_length-ri),,]<-read.nc(f)$GHG_BLOB[li:length(lons),bi:ti,]
# gh_blob_sub[(lon_length-ri+1):lon_length,,]<-read.nc(f)$GHG_BLOB[1:ri,bi:ti,]
# close.nc(f)
# 
# #load the 500 mb Z heights
# z_hgt_sub<-array(NA,c(lon_length,lat_length,time_length))
# tindex=1
# for (x in c("06","07","08")){
#   f<-open.nc(paste(files_dir,"/ERA_2003/ERA_2003_",x,"_vars_z500.nc",sep=""))
#   ntime<-length(read.nc(f)$time)
#   tlast=tindex+(ntime-1)
#   z_hgt_sub[1:(lon_length-ri),,tindex:tlast]<-read.nc(f)$Z[li:length(lons),bi:ti,]
#   z_hgt_sub[(lon_length-ri+1):lon_length,,tindex:tlast]<-read.nc(f)$Z[1:ri,bi:ti,]
#   close.nc(f)
#   tindex=tlast+1
# }
# 
# #7/26 0Z only GHG blob
# #7/31 6Z all 3
# 
# #NOTE: recently downloaded surface data file and corresponding
# # older files have different time units!!!
# time_format<-as.Date(timevar/24, origin="1900-01-01")
# time_hours<-timevar%%24
# first_ind<-which(time_format=="2003-07-01")
# 
# 
# temp.cols<-colorRampPalette(c("purple","blue","cyan4","green","yellow","orange","red", "darkred"))(28)

load("~/block_r_data/eur_heat_wave.RData")
for (t in first_ind:length(time_format)){
#for(t in first_ind:first_ind+3){
  fname<-sprintf("~/Dropbox/eur_heat_wave/tcol/eur_heat_%s_%02dZ.png",time_format[t],time_hours[t])
  png(fname,height=600,width=800)
  filled.contour(lon_seq,lat_seq,t_sub_c[,,t],main=sprintf("2m Temp %s %02dZ, PV (green) Z (blue)",time_format[t],time_hours[t]),
    levels=seq(-65,55,5),col=temp.cols,plot.axes={axis(1);axis(2); 
    map('world',add=TRUE,xlim=c(-110,50),ylim=c(25,75));  
      contour(lon_seq,lat_seq,z_hgt_sub[,,t],levels=seq(4500,6100,50),add=TRUE);
    contour(lon_seq,lat_seq,pv_blob_sub[,,t],levels=c(0,1),add=TRUE,col="chartreuse4",drawlabels=FALSE,lwd=5);
      contour(lon_seq,lat_seq,z_blob_sub[,,t],levels=c(0,1),add=TRUE,col="cornflowerblue",drawlabels=FALSE,lwd=5);
      contour(lon_seq,lat_seq,gh_blob_sub[,,t],levels=c(0,1),add=TRUE,col="purple",drawlabels=FALSE,lwd=5)})  
    dev.off()
}

hgt.cols<-colorRampPalette(c("purple","blue","cyan4","green","yellow", "orange","red","darkred"))(34)
for (t in first_ind:length(time_format)){
#for(t in first_ind:first_ind+3){
  fname<-sprintf("~/Dropbox/eur_heat_wave/noT/eur_heat_%s_%02dZ_noT.png",time_format[t],time_hours[t])
  png(fname,height=600,width=800)

 map('world',xlim=c(-110,50),ylim=c(25,75),fill=TRUE,col="grey")
 title(sprintf("Z500 %s %02dZ, PV* (green) Z* (blue) ZG (purple)",time_format[t],time_hours[t]))
 map.axes()
 contour(lon_seq,lat_seq,z_hgt_sub[,,t],levels=seq(4500,6100,50),drawlabels=FALSE,add=TRUE,col=hgt.cols,lwd=2)
 contour(lon_seq,lat_seq,pv_blob_sub[,,t],levels=c(0,1),add=TRUE,col="chartreuse4",drawlabels=FALSE,lwd=5)
 contour(lon_seq,lat_seq,z_blob_sub[,,t],levels=c(0,1),add=TRUE,col="cornflowerblue",drawlabels=FALSE,lwd=5)
 contour(lon_seq,lat_seq,gh_blob_sub[,,t],levels=c(0,1),add=TRUE,col="purple",drawlabels=FALSE,lwd=5)
  dev.off()
}
#combine 
