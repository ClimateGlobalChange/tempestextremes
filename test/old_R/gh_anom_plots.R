source("~/tempestextremes/test/sector_funcs.R")

#Creating a plot for European heat wave
#Show the surface temperature pattern with the corresponding blocking patterns

#First: anomalies per season per sector
files_dir<-"/Volumes/ExFAT_drive/ERA_files"
img_dir<-"/Users/mariellep/figs"

sectors<-c("NA","NC","NP","SA","SI","SP")
colors<-c("red","orange","pink","blue","green","purple")
seasons<-c("DJF","MAM","JJA","SON")

#datas<-c("climo","2xCO2","SSTplus2")
LEFT_BOUND=c(250, 30, 130, 290, 20, 120)
RIGHT_BOUND=c(50, 150, 270, 40, 140, 310)
MIN_LAT=c(25, 25, 25, -75, -75, -75)
MAX_LAT=c(75, 75, 75, -25, -25, -25)

ystart<-1981
yend<-2000
mstart<-12

#Get array indices for bounds
li<-rep(NA,6)
ri<-rep(NA,6)
bi<-rep(NA,6)
ti<-rep(NA,6)

f<-open.nc("/Volumes/ExFAT_drive/ERA_files/ERA_1980/ERA_1980_01_vars_z500_devs.nc")
lats<-read.nc(f)$lat
lons<-read.nc(f)$lon
close.nc(f)

for (d in 1:6){
  li[d]<-get_index(lons,LEFT_BOUND[d])
  ri[d]<-get_index(lons,RIGHT_BOUND[d])
  bi[d]<-get_index(lats,MIN_LAT[d])
  ti[d]<-get_index(lats,MAX_LAT[d])
}

for (s in 1:length(seasons)){
  #grab each file for season
  if (seasons[s]=="DJF"){
    f1name<-sprintf("%s/ERA_%04d/ERA_%04d_%02d_vars_z500_devs.nc",files_dir,ystart-1,ystart-1,mstart)
    mstart=1
  }else{
    f1name<-sprintf("%s/ERA_%04d/ERA_%04d_%02d_vars_z500_devs.nc",files_dir,ystart,ystart,mstart)
    mstart=mstart+1
  }
  f1<-open.nc(f1name)
  devs1<-read.nc(f1)$ADGH
  close.nc(f1)
  f2name<-sprintf("%s/ERA_%04d/ERA_%04d_%02d_vars_z500_devs.nc",files_dir,ystart,ystart,mstart)
  f2<-open.nc(f2name)
  mstart=mstart+1
  f3name<-sprintf("%s/ERA_%04d/ERA_%04d_%02d_vars_z500_devs.nc",files_dir,ystart,ystart,mstart)
  f3<-open.nc(f3name)
  mstart=mstart+1
  devs2<-read.nc(f2)$ADGH
  close.nc(f2)
  devs3<-read.nc(f3)$ADGH
  close.nc(f3)
  #Do hist per sector
  #dvec<-c()
  for (n in 1:length(sectors)){
    dvec<-c()
    #atlantic sector, grab info 
    if (sectors[n]=="NA" | sectors[n]=="SA"){
      subs1<-devs1[ri[n]:length(lons),ti[n]:bi[n],]
      subs2<-devs2[ri[n]:length(lons),ti[n]:bi[n],]
      subs3<-devs3[ri[n]:length(lons),ti[n]:bi[n],]
      dvec<-c(subs1,subs2,subs3)
      subs1<-devs1[1:li[n],ti[n]:bi[n],]
      subs2<-devs1[1:li[n],ti[n]:bi[n],]
      subs3<-devs1[1:li[n],ti[n]:bi[n],]
      dvec<-c(dvec,subs1,subs2,subs3)
    }else{
      subs1<-devs1[li[n]:ri[n],bi[n]:ti[n],]
      subs2<-devs1[li[n]:ri[n],bi[n]:ti[n],]
      subs3<-devs1[li[n]:ri[n],bi[n]:ti[n],]
      dvec<-c(subs1,subs2,subs3)
    }

    dsub<-dvec[which(dvec>0)]
    dsub_plot<-dsub
    perc<-quantile(dsub,.9)
    perc_plot<-perc


    if (n==1){
      png(sprintf("%s/%s_Zanom_dist.png",img_dir,seasons[s]),width=800,height=600)
      plot(density(dsub_plot),main=seasons[s],col=colors[n],lwd=3)
    }else{
      lines(density(dsub_plot),col=colors[n],lwd=3)
    }
    abline(v=perc_plot,col=colors[n],lwd=2)
    if (n==6){
      abline(v=170,lwd=2,lty=2)
      dev.off()
    }
    mesg<-sprintf("%s %s: top 10 pct is %f",seasons[s],sectors[n],perc)
    print(mesg)
  }
  
}