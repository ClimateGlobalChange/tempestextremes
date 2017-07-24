#Test script to overlap blob densities from sector limited blocking

source("~/tempestextremes/test/sector_funcs.R")

box_out<-function(){
  #Boxes defining sectors
  #NH ATL
  rect(260,25,365,75,lwd=2)
  text(((260+360)/2)+5,50,"NA",cex=2,font=2)
  rect(0,25,40,75,lwd=2)
  text(20,50,"NA",cex=2,font=2)
  #NH PAC
  rect(140,25,260,75,lwd=2)
  text((140+260)/2,50,"NP",cex=2,font=2)
  
  #CONTINENTS
  rect(40,25,140,75,lwd=2)
  text(100,50,"NC",cex=2,font=2)
  
  #SH ATL
  rect(300,-75,365,-25,lwd=2)
  text((300+360)/2,-50,"SA",cex=2,font=2)
  rect(-5,-75,30,-25,lwd=2)
  text(15,-50,"SA",cex=2,font=2)
  
  #SH PAC
  rect(130,-75,300,-25,lwd=2)
  text((130+300)/2,-50,"SP",cex=2,font=2)
  
  #Indian Ocean
  rect(30,-75,130,-25,lwd=2)
  text(80,-50,"SI",cex=2,font=2)
}
sectors<-c("NA","NC","NP","SA","SI","SP")
seasons<-c("DJF","MAM","JJA","SON")
datas<-c("climo","2xCO2","SSTplus2")
LEFT_BOUND=c(250, 30, 130, 290, 20, 120)
RIGHT_BOUND=c(50, 150, 270, 40, 140, 310)
MIN_LAT=c(25, 25, 25, -75, -75, -75)
MAX_LAT=c(75, 75, 75, -25, -25, -25)

d.cols<-colorRampPalette(c("darkgreen", "green","white","magenta", "purple4"))(12)
d.cols[6:7] = "#FFFFFF"


fdir="/Users/mariellep/ERA_CLIMO_DIFF/GH_minus_PV"
setwd(fdir)
  
for (dat in datas){
  patt=paste(dat,"*dens*.nc",sep="")
  print(patt)
  files_list=list.files(pattern=glob2rx(patt))
  for (season in seasons){
    #Subsetting files list by season
    files_sub<-files_list[grep(season,files_list)]
    #Initialize array
    #get lats,lons
    f1<-open.nc(files_sub)
    lats<-read.nc(f1)$lat
    lons<-read.nc(f1)$lon
    d<-read.nc(f1)$dens
    close.nc(f1)
    

    
    png(paste(dat,"_GHmPV_diff_",season,".png",sep=""),width=1000,height=600)
    filled.contour(lons,lats,d,main=paste("GH-PV",dat,season),levels=seq(-.12,.12,.02),col= d.cols,
                   plot.axes={axis(1);axis(2); map('world2',add=TRUE,cex.lab=3); box_out()})
    dev.off()

    

    
    
  }
}
  










