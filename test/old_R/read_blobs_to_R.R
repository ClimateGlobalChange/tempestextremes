library(RNetCDF)

source("~/tempestextremes/test/sector_funcs.R")
fdir<-"/Volumes/ExFAT_drive/ERA_files/"
sectors<-c("NA","NC","NP","SA","SI","SP")
LEFT_BOUND=c(250, 30, 130, 290, 20, 120)
RIGHT_BOUND=c(50, 150, 270, 40, 140, 310)
MIN_LAT=c(25, 25, 25, -75, -75, -75)
MAX_LAT=c(75, 75, 75, -25, -25, -25)

#Axes
f<-open.nc("/Volumes/ExFAT_drive/ERA_files/ERA_avg/ERA_1980_2005_pv_avgDFT.nc")
lats<-read.nc(f)$lat
lons<-read.nc(f)$lon
close.nc(f)

lons_c<-ifelse(lons<181,lons,lons-360)

subdir<-"ERA_blobs/"
for (i in 1:6){
  seci<-which(sectors==sectors[i])
  lb<-get_index(lons,LEFT_BOUND[seci])
  rb<-get_index(lons,RIGHT_BOUND[seci])
  if (lb>rb){
    lons_sub<-c(seq(lb,length(lons)),seq(1,rb))
  }else{
    lons_sub<-seq(lb,rb)
  }
  
  tb<-get_index(lats,MIN_LAT[seci])
  bb<-get_index(lats,MAX_LAT[seci])
  if (tb>bb){
    lats_sub<-seq(bb,tb)
  }else{
    lats_sub<-seq(tb,bb)
  }
  
  latsize<-length(lats_sub)
  lonsize<-length(lons_sub)
  
  
  #Size of time axis:
  tsz<-4*365*26
  
  pv_anom<-array(NA,c(lonsize,latsize,tsz))
  z_anom<-array(NA,c(lonsize,latsize,tsz))
  ghg<-array(NA,c(lonsize,latsize,tsz))
  hrs<-c()
  start_t<-1

  #READ IN ALL DATA FOR WORKSPACE
  for (y in 1980:2005){

    for (m in c("MAM","JJA","SON","DJF")){
      fname_pvanom<-paste(fdir,subdir,sprintf("ERA_%04d_%s_%s_blobs.nc",y,m,sectors[i]),sep="")
      #READ PV ANOM DATA
      print(sprintf("Reading %s",fname_pvanom))
      #Open PV
      fpv<-open.nc(fname_pvanom)
      #Add hours to time axis
      t<-read.nc(fpv)$time
      tlen<-length(t)
      hrs<-c(hrs,t)
      #Get end index for time axis
      end_t<-start_t+tlen-1
      #Save PV data
      pv_anom[,,start_t:end_t]<-read.nc(fpv)$PV_BLOB[lons_sub,lats_sub,]
      close.nc(fpv)
      #READ Z ANOM DATA
      fname_zanom<-paste(fdir,subdir,sprintf("ERA_%04d_%s_%s_Zblobs.nc",y,m,sectors[i]),sep="")
      print(sprintf("Reading %s",fname_zanom))
      #Open Z anom
      fz<-open.nc(fname_zanom)
      z_anom[,,start_t:end_t]<-read.nc(fz)$Z_BLOB[lons_sub,lats_sub,]
      close.nc(fz)
      
      #READ Z ANOM DATA
      fname_ghg<-paste(fdir,subdir,sprintf("ERA_%04d_%s_%s_GHGblobs.nc",y,m,sectors[i]),sep="")
      print(sprintf("Reading %s",fname_ghg))
      #Open Z anom
      fg<-open.nc(fname_ghg)
      ghg[,,start_t:end_t]<-read.nc(fg)$GHG_BLOB[lons_sub,lats_sub,1:tlen]
      close.nc(fg)
      #Increment next t
      start_t<-end_t+1
    }
  }
  image_name<-paste(sprintf("~/block_r_data/%s_pv_z_ghg_block_data.RData",sectors[i]))
  time_format<-as.Date(hrs/24, origin="1800-01-01")
  time_hrs<-hrs%%24
  lons_seq<-lons[lons_sub]
  lons_seq_c<-lons_c[lons_sub]
  lats_seq<-lats[lats_sub]
  save(list=c("pv_anom","z_anom","hrs","lats_sub","ghg","lats_seq","lons_seq","lons_seq_c",
              "lons_sub","lats","lons","time_format","time_hrs"),file=image_name)
}


