library("RNetCDF")
library("maps")
library("maptools")
library("mapproj")
library("abind")
library("RColorBrewer")
library("MASS")











img_dir=paste("~/data_netcdf/",varname,sep="")
if (!dir.exists(img_dir)){
  dir.create(img_dir)
}

setwd(img_dir)

#################################


nYearCompare=21
nPerCompare=4

data_vec<-c("ERA", "climo", "2xCO2","SSTplus2","SSTplus2_2xCO2")
labs<-c("ERA","climo","2xCO2","SST+2","SST+2_2xCO2")
seasons_vec<-c("DJF","JJA","MAM","SON")
sec=c("ATL","PAC")
hemi=c("NH","SH")

breaks_NH=seq(24,82,1)
breaks_SH=seq(-82,-24,1)
###################################################

# Correlation data: Generates Spearman correlation coefficients
#Correlation matrix between ERA and other datasets
for (season in seasons_vec){
  #For each hemisphere
  for (h in hemi){
    #Calculate breaks for centerlat
    if (varname=="centlat"){
      if (h=="NH"){
        breaks=breaks_NH
      }else{
        breaks=breaks_SH
      }
    }
    for(s in sec){
      print(paste("Spearman correlation between datasets for",varname,season,h,s))
      dat=df_tot[df_tot$Sector==s & df_tot$Hemi==h & df_tot$Season==season,c("Dataset","nYears","centlon_c",varname)]
      #Calculate breaks for centlon
      if (varname=="centlon"){
        if (s=="PAC"){
          varname_plot="centlon"
        }else if (s=="ATL"){
          varname_plot="centlon_c"
        }
      }else if (varname != "centlon"){
        varname_plot=varname
      }
      lower=floor(min(dat[,varname_plot]))
      upper=ceiling(max(dat[,varname_plot]))
      break_interval=(upper-lower)/60
      breaks=seq(lower,upper,break_interval)
      
      mat=data.frame(matrix(nrow=length(breaks)-1,ncol=length(data_vec)))
      colnames(mat)<-data_vec
      
      for (d in 1:length(data_vec)){
        sub=dat[dat$Dataset==data_vec[d],]
        if (data_vec[d]=="ERA"){
          t=4
        }else{
          t=8
        }
        
        n=sub$nYears[1]
        scale=(nYearCompare*nPerCompare)/(t*n)
        h1=hist(sub[,varname_plot],breaks=breaks,plot=FALSE,include.lowest=TRUE)
        mat[,d]=h1$counts*scale
        
      }
      print(cor(mat,method="pearson"))
    }
  }
}

###################################################

#Outputs Q-Q plots (climo vs others)
for (season in seasons_vec){
  for (h in hemi){
    if (varname =="centlat"){
      if (h=="NH"){
        breaks=breaks_NH
      }else{
        breaks=breaks_SH
      }
    }
    for(s in sec){
      dat=df_tot[df_tot$Sector==s & df_tot$Hemi==h & df_tot$Season==season,c("centlon_c","Dataset","nYears",varname)]
      #Get the extent of the center longitude values 
      if (varname=="centlon"){
        if (s=="PAC"){
          varname_plot="centlon"
        }else if (s=="ATL"){
          varname_plot="centlon_c"
        }
      }else if (varname != "centlon"){
        varname_plot=varname
      }
      lower=floor(min(dat[,varname_plot]))
      upper=ceiling(max(dat[,varname_plot]))
      break_interval=(upper-lower)/60
      breaks=seq(lower,upper,break_interval)
      
      mat=data.frame(matrix(nrow=length(breaks)-1,ncol=length(data_vec)))
      colnames(mat)<-data_vec
      for (d in 1:length(data_vec)){
        sub=dat[dat$Dataset==data_vec[d],]
        if (data_vec[d]=="ERA"){
          t=4
        }else{
          t=8
        }
        n=sub$nYears[1]
        scale=(nYearCompare*nPerCompare)/(t*n)
        h1=hist(sub[,varname_plot],breaks=breaks,plot=FALSE,include.lowest=TRUE)
        mat[,d]=h1$counts*scale
        
      }
      pname=paste(season,"_",h,"_",s,"_",varname,"_climo_compare.png",sep="")
      png(pname,width=1200,height=300)
      par(mfrow=c(1,4))
      par(mar=c(4,4,1,1)+.1)
      for (i in 2:2){
        for (j in 1:length(data_vec)){
          if (i!=j){
            qqplot(mat[,j],mat[,i],xlab=labs[j],ylab=labs[i],cex.lab=1.6)
            abline(a=0,b=1)
          }
        }
      }
      dev.off()
    }
    
  }
}

###################################################

#Outputs histograms of the data
#Loop through datasets and seasons and plot histogram of center latitudes

#Histograms need to be "adjusted" because of different counts due to diff # years and time steps!!!

for (season in seasons_vec){
  for (h in hemi){
    if (varname=="centlat"){
      if (h=="NH"){
        breaks=breaks_NH
      }else{
        breaks=breaks_SH
      }
    }
    for (s in sec){
      dat=df_tot[df_tot$Sector==s & df_tot$Hemi==h & df_tot$Season==season,c("centlon_c","Dataset","nYears",varname)]
      if (varname=="centlon"){
        if (s=="PAC"){
          varname_plot="centlon"
        }else if (s=="ATL"){
          varname_plot="centlon_c"
        }
      }else if (varname != "centlon"){
        varname_plot=varname
      }
      
      lower=min(dat[,varname_plot])
      upper=max(dat[,varname_plot])
      break_interval=(upper-lower)/60
      
      breaks=seq(lower,upper,break_interval)
      
      for (dataname in data_vec){
        #Define histogram color
        if (dataname=="climo"){
          histcol="blue"
        }else if (dataname=="2xCO2"){
          histcol="cyan4"
        }else if (dataname=="SSTplus2"){
          histcol="blueviolet"
        }else if (dataname=="SSTplus2_2xCO2"){
          histcol="darkblue"
        }else if (dataname=="ERA"){
          histcol="aquamarine2"
        }
        #Number of time steps
        if(dataname == "ERA"){
          nPer=4
        }else{
          nPer=8
        }
        histname=paste(season,"_",h,"_",s,"_",varname,"_hist_",dataname,".png",sep="")
        png(filename=histname,width=600,height=300)
        xlim=range(breaks)
        mt=paste(dataname,"blocking",varname,"counts for",season,h,s)
        sub=dat[dat$Dataset==dataname,]
        
        nYears=sub$nYears[1]
        
        mult_factor = (nYearCompare*nPerCompare)/(nYears*nPer)
        
        hist_data=hist(sub[,varname_plot],breaks=breaks,plot=FALSE)
        hist_data$counts=round(hist_data$counts*mult_factor)
        plot(hist_data, main=mt,xlab=varname,col=histcol,xlim=xlim,ylim=ylim,cex.main=1.3)
        dev.off()
      }
    }
  }
}

#   #m<-map("world",fill=TRUE,col="grey",mar=rep(.5,4))
#   m<-map("world",mar=rep(.5,4))
#   k <- 9
#   my.cols <- rev(brewer.pal(k, "YlGn"))
#   z <- kde2d(df$centlon_c, df$centlat, n=50)
#   #points(df$centlon_c,df$centlat, xlab="lon", ylab="lat", pch=19, cex=.4)
#   contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE,lwd=2)
#   rect(140,25,180,75,border="red")
#   rect(-180,25,-100,75,border="red")
#   #rect(40,25,140,75,border="blue")
#   rect(-80,25,40,75,border="red")
#   rect(-60,-75,30,-25,border="red")
#   rect(130,-75,180,-25,border="red")
#   rect(-180,-75,-60,-25,border="red")





#################
#NETCDF STUFF
#files go from 02-12-01 to 24-02-25 (ends on 02-28)
listfiles=readLines("test_climo.txt")
nFiles=length(listfiles)

ref=open.nc(listfiles[1])
lat=var.get.nc(ref,"lat")
lon=var.get.nc(ref,"lon")
tunit=att.get.nc(ref,"time","units")
close.nc(ref)

latval=52.3
lonval=198.75
a=get_index(lat,latval)
b=get_index(lon,lonval)

dat=c()
tdf=data.frame()
#Test plots: plot single value
nPerYear=365*8
for (l in listfiles){
  print(l)
  f=open.nc(l)
  t=var.get.nc(f,"time")
  tconv=utcal.nc(tunit,t)
  tlen=length(t)
  v=var.get.nc(f,"INT_ADIPV",start=c(b,a,NA),count=c(1,1,NA))
  close.nc(f)
  dat=append(dat,v)
  tdf=rbind(tdf,tconv)
}

#DJF
#Dates
tdjf=tdf[((tdf$month==12)||(tdf$month==1)||(tdf$month==2)),]
indices=as.numeric(rownames(tdjf))
#Anomaly values (integer value)
djfdat=dat[indices]
days=day_in_season(tdjf)

plot(days,djfdat)



scatterhist = function(x, y, xlab="", ylab=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  plot(x,y)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
        at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
        at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}


#FOR REFERENCE: BINNING
#Define the lats for the longitude bands
#   lat_breaks=seq(floor(min(sub$centlat))-0.5,ceiling(max(sub$centlat))+0.5)
#   lats=seq(floor(min(sub$centlat)),ceiling(max(sub$centlat)))


#   if (sec=="PAC"){
#     lon_vals=sub$centlon_c
#   }
#   else{
#     lon_vals=sub$centlon
#   }
#   lon_breaks=seq(floor(min(lon_vals))-0.5,ceiling(max(lon_vals))+0.5)
#   lons=seq(floor(min(lon_vals)),ceiling(max(lon_vals)))




# 
#   for (i in 1:5){
#     subdf=sub[sub$time==i,c("centlat","centlon")]
#     #bin the lon values
#     lon_bins=as.numeric(cut(subdf$centlon,lon_breaks,labels=FALSE))
#     lon_uniq=unique(lon_bins)
#     for (b in lon_uniq){
#       match=lon_bins[lon_bins==b]
#       if (length(match)>1){
#         #The values need to be averaged
#       }else{
#         
#       }
#     }
#     #bin the lat values and average the lon values
# #     lat_bins=cut(subdf$centlat,lat_breaks,labels=FALSE)
# #     lat_ind=as.numeric(lat_bins)
#     
#   }

#Plot mean lat over time per basin
meanlat_over_time<-function(hemi,sec,col){
  sub=df[(df$Sector==sec & df$Hemi==hemi),]
  t=seq(min(sub$time),max(sub$time))
  tday=t/8+1
  meanlats<-c()
  for (i in t){
    subdf=sub[sub$time==i,"centlat"]
    meanlats<-c(meanlats,mean(subdf))
  }
  #header=paste("Mean center latitude in",hemi,sec,"basin")
  #plot(t,meanlats,type="l",main=header,xlab="time",ylab="lat",xlim=xlim,ylim=ylim)
  lines(tday,meanlats,type="l",col=col)
  
}
trange=range(df$time/8 + 1)

plot(0,type="n",xlim=trange, ylim=c(48,62),
     main="Mean center latitude in NH for DJF",xlab="day in season",ylab="center lat")
meanlat_over_time("NH","PAC","purple")
meanlat_over_time("NH","ATL","blue")
dev.print(device = png, filename = 'DJF_NH_centlat_time.png', width = 1000, height = 500) 

plot(0,type="n",xlim=trange, ylim=c(-70,-54),
     main="Mean center latitude in SH for DJF",xlab="day in season",ylab="center lat")
meanlat_over_time("SH","PAC","purple")
meanlat_over_time("SH","ATL","blue")
dev.print(device = png, filename = 'DJF_SH_centlat_time.png', width = 1000, height = 500) 

#Today's goal: Zonally averaged histogram of center lon values
#i.e. for each longitude band, what is the average center longitude?
zon_avg_hist<-function(hemi,sec){
  sub=df[(df$Sector==sec & df$Hemi==hemi),]
  #Create new lat column
  sub$lat_int=floor(sub$centlat)
  #Now, for each time step, what are the repeat lats?
  t=seq(min(sub$time),max(sub$time))
  centlon_hist_vals<-c()
  for (i in t){
    subdf=sub[sub$time==i,]
    lat_uniq=unique(subdf$lat_int)
    for (l in lat_uniq){
      lon_vals=subdf[subdf$lat_int==l,"centlon"]
      centlon_hist_vals=c(centlon_hist_vals,mean(lon_vals))
    }
  }
  
  
}

#smoothScatter(df$centlon_c,df$centlat)
#Block centerpoints by sector
# abline(v=c(40,140,260,280),col="red")
# text(20,2000,"ATL")
# text(200,2000,"PAC")
# text(270,2000,"CONT1")
# text(310,2000,"ATL")
# text(90,2000,"CONT2")