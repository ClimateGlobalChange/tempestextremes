library(RNetCDF)

source("~/tempestextremes/test/sector_funcs.R")
args=commandArgs(trailingOnly=TRUE)
if (length(args)<3){
   stop("Argument 1: working directory Argument 2: start year Argument 3: end year",call.=FALSE)
}

fdir<-args[1]
setwd(fdir)
data="MERRA"
syear=as.numeric(args[2])
eyear=as.numeric(args[3])

sectors<-c("NA","NC","NP","SA","SI","SP")
if (data=="ERA"){
LEFT_BOUND=c(250, 30, 130, 290, 20, 120)
RIGHT_BOUND=c(50, 150, 270, 40, 140, 310)
MIN_LAT=c(25, 25, 25, -75, -75, -75)
MAX_LAT=c(75, 75, 75, -25, -25, -25)
}

if (data=="MERRA"){
LEFT_BOUND=c(-110, 30, 130, -70, 20, 120)
RIGHT_BOUND=c(50, 150, -90, 40, 140, -50)
MIN_LAT=c(25, 25, 25, -75, -75, -75)
MAX_LAT=c(75, 75, 75, -25, -25, -25)

}

#Axes
f<-open.nc(sprintf("%s/%s_1980_2005_Z500_avg_DFT.nc",fdir,data))
lats<-read.nc(f)$lat
lons<-read.nc(f)$lon
close.nc(f)

nyears=(eyear-syear)+1
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
  tsz<-4*365*nyears
  
#  pv_inst<-array(NA,c(lonsize,latsize,tsz))
  z_inst<-array(NA,c(lonsize,latsize,tsz))
  
  time_format<-c()
  hrs<-c()
  start_t<-1
  #READ IN ALL DATA FOR WORKSPACE
  for (y in syear:eyear){
    subdir<-sprintf("%s_%04d",data,y)
    for (m in 1:12){
      #List the inst files
      plist=list.files(path=subdir,full.names=TRUE,pattern=glob2rx(sprintf("*%04d_%02d*integ.nc",y,m)))
      zlist=list.files(path=subdir,full.names=TRUE,pattern=glob2rx(sprintf("*%04d_%02d*vars.nc4",y,m)))
      nfiles=length(plist) 
      for (d in 1:nfiles){
        #READ Z ANOM DATA
        # fname_zanom<-paste(fdir,subdir,sprintf("%s_%04d_%02d_vars_z500_devs.nc",data,y,m),sep="")
        # fname_pvanom<-paste(fdir,subdir,sprintf("%s_%04d_%02d_vars_integ_devs.nc",data,y,m),sep="")
        fname_zanom<-zlist[d]
 #       fname_pvanom<-plist[d]
        print(sprintf("Reading %s",fname_zanom))
        #Open Z anom
        fz<-open.nc(fname_zanom)
        fp<-read.nc(fz)$lev
        pindex<-which(fp==500)
        #READ PV ANOM DATA
#        print(sprintf("Reading %s",fname_pvanom))
        #Open PV
#        fpv<-open.nc(fname_pvanom)
        #Add hours to time axis
        t<-read.nc(fz)$time
        hrs<-c(hrs,t)
        time_date<-as.Date(t/(24*60), origin=sprintf("%04d-%02d-%02d",y,m,d))
        time_hours<-(t/60)%%24
        time_string<-sprintf("%s_%02d",time_date,time_hours)
        time_format<-c(time_format,time_string)
        #Get end index for time axis
        end_t<-start_t+length(t)-1
        #Save PV data
 #       pv_inst[,,start_t:end_t]<-read.nc(fpv)$VPV[lons_sub,lats_sub,]
 #       close.nc(fpv)
     #   print(sprintf("Reading %s",fname_zanom))
        #Open Z anom
     #   fz<-open.nc(fname_zanom)
     #   fp<-read.nc(fz)$lev
     #   pindex<-which(fp==500)
        z_inst[,,start_t:end_t]<-read.nc(fz)$Z[lons_sub,lats_sub,pindex,]
        close.nc(fz)
        #Increment next t
        start_t<-end_t+1
      }
    }
  }
  image_name<-paste(sprintf("~/block_r_data/%s_%s_z_inst_data.RData",data,sectors[i]))
  lons_seq<-lons[lons_sub]
  lats_seq<-lats[lats_sub]
  save(list=c("z_inst","hrs","lats_sub","lats_seq","lons_seq",
              "lons_sub","lats","lons","time_format"),file=image_name)
}


