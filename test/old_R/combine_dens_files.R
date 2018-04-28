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

blob.cols<-colorRampPalette(c("white","mediumpurple4","blue"))(9)

for (dir in c("GH")){
  fdir=paste("/Users/mariellep/ERA_CLIMO_DIFF/",dir,sep="")
  setwd(fdir)
  
  for (dat in datas){
    patt=paste(dat,"*dens*.nc",sep="")
    print(patt)
    fpath=paste(fdir,"sector_dens",sep="/")
    print(fpath)
    files_list=list.files(pattern=glob2rx(patt),path=fpath)
    setwd(fpath)
    for (season in seasons){
      #Subsetting files list by season
      files_sub<-files_list[grep(season,files_list)]
      #Initialize array
      #get lats,lons
      f1<-open.nc(files_list[1])
      lats<-read.nc(f1)$lat
      lons<-read.nc(f1)$lon
      
      #open new file for writing
      new_fname=paste(fdir,"/",dat,"_",season,"_avg_zdens_time.nc",sep="")
      print(new_fname)
      f2<-create.nc(new_fname)
      dim.def.nc(f2,"lat",length(lats))
      dim.def.nc(f2,"lon",length(lons))

      var.def.nc(f2,"lat","NC_DOUBLE","lat")
      var.put.nc(f2,"lat",lats,1,length(lats))
      att.copy.nc(f1,"lat","long_name",f2,"lat")
      att.copy.nc(f1,"lat","units",f2,"lat")


      var.def.nc(f2,"lon","NC_DOUBLE","lon")
      var.put.nc(f2,"lon",lons,1,length(lons))
      att.copy.nc(f1,"lon","long_name",f2,"lon")
      att.copy.nc(f1,"lon","units",f2,"lon")
      
      close.nc(f1)
      
      #Create an empty array
      ar<-array(rep(NaN, length(sectors)*length(lats)*length(lons)),
                dim=c(length(sectors),length(lons),length(lats)))
      
      d=1
      li<-rep(NA,6)
      ri<-rep(NA,6)
      bi<-rep(NA,6)
      ti<-rep(NA,6)
      for (s in sectors){
        #File name
        fname=files_sub[grep(s,files_sub)]
        #open file and read values into array
        nc<-open.nc(fname)
        dens<-read.nc(nc)$dens
        ar[d,,]=dens
        
        #Indices for bounds
        li[d]<-get_index(lons,LEFT_BOUND[d])
        ri[d]<-get_index(lons,RIGHT_BOUND[d])
        bi[d]<-get_index(lats,MIN_LAT[d])
        ti[d]<-get_index(lats,MAX_LAT[d])
        d=d+1
        close.nc(nc)
      }
      
      #Need to look at the overlaps
      ar2<-array(rep(0,length(lons)*length(lats)),dim=c(length(lons),length(lats)))
      #Start filling in densities
      #NA; account for periodic condition
      ar2[ri[3]:length(lons),bi[1]:ti[1]]=ar[1,ri[3]:length(lons),bi[1]:ti[1]]
      ar2[1:li[2],bi[1]:ti[1]]=ar[1,1:li[2],bi[1]:ti[1]]
      #NC
      ar2[ri[1]:li[3],bi[1]:ti[1]]=ar[2,ri[1]:li[3],bi[1]:ti[3]]
      #NP
      ar2[ri[2]:li[1],bi[1]:ti[1]]=ar[3,ri[2]:li[1],bi[1]:ti[1]]
      #SA
      ar2[ri[6]:length(lons),bi[4]:ti[4]]=ar[4,ri[6]:length(lons),bi[4]:ti[4]]
      ar2[1:li[5],bi[4]:ti[4]]=ar[4,1:li[5],bi[4]:ti[4]]
      #SI
      ar2[ri[4]:li[6],bi[4]:ti[4]]=ar[5,ri[4]:li[6],bi[4]:ti[4]]
      #SP
      ar2[ri[5]:li[4],bi[4]:ti[4]]=ar[6,ri[5]:li[4],bi[4]:ti[4]]
      
      ###############
      #Pairs: 
      #NC L and NA R
      N1_ar<-ar[1,li[2]:ri[1],bi[1]:ti[1]]
      N2_ar<-ar[2,li[2]:ri[1],bi[1]:ti[1]]
      NC_NA=(N1_ar + N2_ar)/2
      ar2[li[2]:ri[1],bi[1]:ti[1]]=NC_NA
      #cat("NC:",N2_ar[17,52],",NA:",N1_ar[17,52],"avg:",NC_NA[17,52])
      #NP L, NC R
      N1_ar<-ar[2,li[3]:ri[2],bi[1]:ti[1]]
      N2_ar<-ar[3,li[3]:ri[2],bi[1]:ti[1]]
      NP_NC=(N1_ar + N2_ar)/2
      ar2[li[3]:ri[2],bi[1]:ti[1]]=NP_NC
      #NA L, NP R
      N1_ar<-ar[1,li[1]:ri[3],bi[1]:ti[1]]
      N2_ar<-ar[3,li[1]:ri[3],bi[1]:ti[1]]
      NA_NP=(N1_ar + N2_ar)/2
      ar2[li[1]:ri[3],bi[1]:ti[1]]=NA_NP
      #SI L, SA R
      N1_ar<-ar[5,li[5]:ri[4],bi[4]:ti[4]]
      N2_ar<-ar[4,li[5]:ri[4],bi[4]:ti[4]]
      SI_SA=(N1_ar + N2_ar)/2
      ar2[li[5]:ri[4],bi[4]:ti[4]]=SI_SA
      #cat("SA:",N2_ar[1,9],",SI:",N1_ar[1,9],"avg:",SI_SA[1,9])
      
      #SP L, SI R
      N1_ar<-ar[6,li[6]:ri[5],bi[4]:ti[4]]
      N2_ar<-ar[5,li[6]:ri[5],bi[4]:ti[4]]
      SP_SI=(N1_ar + N2_ar)/2
      ar2[li[6]:ri[5],bi[4]:ti[4]]=SP_SI
      #cat("SI:",N2_ar[1,9],",SP:",N1_ar[1,9],"avg:",SP_SI[1,9])
      
      #SA L, SP R
      N1_ar<-ar[4,li[4]:ri[6],bi[4]:ti[4]]
      N2_ar<-ar[6,li[4]:ri[6],bi[4]:ti[4]]
      SA_SP=(N1_ar + N2_ar)/2
      ar2[li[4]:ri[6],bi[4]:ti[4]]=SA_SP
      
      # png(paste(dat,"_",dir,"_dens_abs_",season,".png",sep=""),width=1000,height=600)
      # filled.contour(lons,lats,ar2,main=paste(dir,dat,season),levels=seq(0,.18,.02),col= blob.cols,
      #                plot.axes={axis(1);axis(2); map('world2',add=TRUE,cex.lab=3); box_out()})
      # dev.off()
      #cat("SP:",N2_ar[1,9],",SA:",N1_ar[1,9],"avg:",SA_SP[1,9])
      
      #Now save this as a new netCDF!!
      var.def.nc(f2,"dens","NC_DOUBLE",c("lon","lat"))
      var.put.nc(f2,"dens",ar2,c(1,1),c(length(lons),length(lats)))
      close.nc(f2)
      
      
    }
  }
  
}









