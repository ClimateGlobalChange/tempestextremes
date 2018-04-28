library(RNetCDF)

source("~/tempestextremes/test/sector_funcs.R")
fdir<-"/Volumes/ExFAT_drive/ERA_files/"
sectors<-c("NA","NC","NP","SA","SI","SP")
LEFT_BOUND=c(250, 30, 130, 290, 20, 120)
RIGHT_BOUND=c(50, 150, 270, 40, 140, 310)
MIN_LAT=c(25, 25, 25, -75, -75, -75)
MAX_LAT=c(75, 75, 75, -25, -25, -25)

#Axes
f<-open.nc("/Volumes/ExFAT_drive/ERA_files/ERA_avg/ERA_1980_2005_z500_avgDFT.nc")
lats<-read.nc(f)$lat
lons<-read.nc(f)$lon
close.nc(f)

for (i in c(2)){
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
  
  pv_inst<-array(NA,c(lonsize,latsize,tsz))
  z_inst<-array(NA,c(lonsize,latsize,tsz))
  
  hrs2<-c()
  start_t<-1
  #READ IN ALL DATA FOR WORKSPACE
  for (y in 1980:2005){
    subdir<-sprintf("ERA_%04d/",y)
    for (m in 1:12){
      fname_pvanom<-paste(fdir,subdir,sprintf("ERA_%04d_%02d_vars_integ.nc",y,m),sep="")
      #READ PV ANOM DATA
      print(sprintf("Reading %s",fname_pvanom))
      #Open PV
      fpv<-open.nc(fname_pvanom)
      #Add hours to time axis
      t<-read.nc(fpv)$time

      #Get end index for time axis
      end_t<-start_t+length(t)-1
      if (m==2){
        end_t<-end_t-4
        hrs2<-c(hrs2,t[1:(length(t)-4)])
        pv_inst[,,start_t:end_t]<-read.nc(fpv)$IPV[lons_sub,lats_sub,1:(length(t)-4)]
      }else{
        hrs2<-c(hrs2,t)
        pv_inst[,,start_t:end_t]<-read.nc(fpv)$IPV[lons_sub,lats_sub,]       
      }
      close.nc(fpv)
      #READ Z ANOM DATA
      fname_zanom<-paste(fdir,subdir,sprintf("ERA_%04d_%02d_vars_z500.nc",y,m),sep="")
      print(sprintf("Reading %s",fname_zanom))
      #Open Z anom
      fz<-open.nc(fname_zanom)
      if (m==2){
        z_inst[,,start_t:end_t]<-read.nc(fz)$Z[lons_sub,lats_sub,1:(length(t)-4)]
      }else{
        z_inst[,,start_t:end_t]<-read.nc(fz)$Z[lons_sub,lats_sub,]
      }

      close.nc(fz)
      #Increment next t
      start_t<-end_t+1
    }
  }
  image_name<-paste(sprintf("~/block_r_data/%s_pv_z_inst_data.RData",sectors[i]))
  time_format2<-as.Date(hrs2/24, origin="1800-01-01")
  time_hrs2<-hrs2%%24
  save(list=c("pv_inst","z_inst","hrs2","lats_sub",
              "lons_sub","lats","lons","time_format2","time_hrs2"),file=image_name)
}


