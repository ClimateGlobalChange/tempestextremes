library("RNetCDF")
library("maps")
library("maptools")
library("mapproj")
library("abind")
library("RColorBrewer")
library("MASS")
setwd("~/data_netcdf")

#Data we have:
## Original variables
## Blobs
## Density of blobs
## Stats: max/min/center lat/lon, area

#NOTE!!!!!!!!!!!!!!!!!!!
#NETCDF VARIABLES ARE INDEXED AS LON,LAT,TIME!!

#Function that returns the index that most closely matches the specified degree
#Arguments: lat or lon array, degree
get_index<-function(arr,n1){
  return(which.min(abs(arr-n1)))
}

#Function to output an array within the bounds of the specified sector
#Arguments: variable, lat array, lon array, sector (string), hemisphere (string)
get_arr<-function(var,lat,lon,sec,hemi){

  if (hemi=="NH"){
    l1=get_index(lat,15.)
    l2=get_index(lat,75.)

    if (sec=="ATL"){
      n1=get_index(lon,280.)
      n2=get_index(lon,40.) 
    }else if (sec=="PAC"){
      n1=get_index(lon,140.)
      n2=get_index(lon,260.)
    }else if (sec=="CONT1"){
      n1=get_index(lon,261.) 
      n2=get_index(lon,279.)
    }else if (sec=="CONT2"){
      n1=get_index(lon,41.)
      n2=get_index(lon,139.)
    }
    print(n1)
    print(n2)
  }else if (hemi=="SH"){
    l1=get_index(lat,-75.)
    l2=get_index(lat,-15.)
    
    if (sec=="ATL"){
      n1=get_index(lon,300.)
      n2=get_index(lon,30.)
    }else if (sec=="PAC"){
      n1=get_index(lon,130.)
      n2=get_index(lon,300.)
    }else if (sec=="IO"){
      n1=get_index(lon,31.)
      n2=get_index(lon,129.)
    }
    
  }
  if (sec=="ATL"){
    ATL1=var[n1:length(lon),l1:l2,]
    ATL2=var[1:n2,l1:l2,]
    outarr=abind(ATL1,ATL2,along=1)
  }else{
    outarr=var[n1:n2,l1:l2,]
  }
  return(outarr)
}

#Function that returns the day # in the season
#Arguments: time array (which contains year, month, day, hour)
day_in_season<-function(tarr){
  ifelse(
    (tarr$month==12) | (tarr$month==3) |(tarr$month==6) | (tarr$month==9),
    tarr$day, 
    ifelse((tarr$month==1) | (tarr$month==4), 
           31+tarr$day, 
           ifelse((tarr$month==7) | (tarr$month==10),
                  30+tarr$day,
                  ifelse(tarr$month==2, 
                         62+tarr$day,
                         61+tarr$day)
           )
    )
  )+tarr$hour/24
}
lon_convert<-function(lon){
  distFrom180=lon-180.
  return(ifelse(
    distFrom180<0,
    lon,
    -(180-distFrom180)
  ))
}

#Today's goal: create a histogram of block centerpoints
#  by basin and season
#Note that this will need to change for ERA...

sector_df<-function(df){
  ifelse(
    df$Hemi=="NH",
    ifelse(
      #ATLANTIC
      ((df$centlon < 40.) | (df$centlon > 280.)),
      "ATL",
      #PACIFIC
      ifelse(
        ((df$centlon > 140.) & (df$centlon < 260.)),
        "PAC",
        ifelse(
          ((df$centlon >=260.) & (df$centlon)<=280.),
          "CONT1",
          "CONT2"
        )
      )
    ),#SH SECTORS
    ifelse(
      #ATLANTIC
      ((df$centlon < 30.) | (df$centlon > 300.)),
      "ATL",
      ifelse(
        #PACIFIC
        ((df$centlon > 130.) & (df$centlon < 300.)),
        "PAC",
        #INDIAN OCEAN
        "IO"
      )
    )
  )
}

listfiles=readLines("DJF_climo_statlist")
nFiles=length(listfiles)

#Find out total number of lines, + header lines
tot_flines=0
tot_hlines=0
for (x in 1:nFiles){
  f=listfiles[x]
  lines=readLines(f)
  nlines=length(lines)
  
  tot_flines = tot_flines + nlines
  #Grab all instances of header line
  BlobCounts=0
  for (i in 1:nlines){
    if (!is.na(pmatch("Blob",lines[i]))){
      BlobCounts = BlobCounts + 1
    }
  }
  tot_hlines = tot_hlines + BlobCounts
}
nrows=tot_flines-tot_hlines

#Initialize data frame
col = c('Season','year','Blob','time','minlat','maxlat','minlon','maxlon','centlat','centlon','area')
ncols=length(col)
df=data.frame(matrix(ncol=ncols, nrow=nrows))
colnames(df)<-col

start_df_index=1
end_df_index=0

#Read in data!
#for(x in 1:1){
for (x in 1:nFiles){
  #Read file
  f=listfiles[x]
  #Get number of lines
  lines=readLines(f)
  nlines=length(lines)
  #Find all instances of header lines
  BlobIndices=c()
  for (i in 1:nlines){
    if (!is.na(pmatch("Blob",lines[i]))){
      BlobIndices=append(BlobIndices,i)
    }
  }
  
  #Grab year and season from file name
  fstring=unlist(strsplit(f,split="/"))
  ystring=unlist(strsplit(fstring[length(fstring)],split="_"))
  season=ystring[1]
  year=as.numeric(ystring[2])

  #split strings and convert to data frame
#  for (n in 1:2){
  for (n in 1:length(BlobIndices)){
    print(n)
    if (n<(length(BlobIndices))){
      i1=BlobIndices[n]
      i2=BlobIndices[n+1]-1
    }else{
      i1=BlobIndices[n]
      i2=length(lines)
    }
    sub=lines[i1:i2]
    #number of lines (minus header)
    lsub=length(sub)-1
    #Number of the blob
    title=unlist(strsplit(sub[1],split=" "))
    nBlob=as.numeric(title[2])

    #Vector with Blob number, year number, season
    bvec=rep(nBlob,lsub)
    yvec=rep(year,lsub)
    svec=rep(season,lsub)
    
    end_df_index=end_df_index + lsub

    #Data
    dat=do.call(rbind,strsplit(sub[2:length(sub)],split="\t"))

    final_dat=cbind(svec,yvec,bvec,dat)

    df[start_df_index:end_df_index,] = final_dat
    start_df_index = end_df_index + 1
  }
}
df[,2:11]=lapply(df[2:11],as.numeric)
df$Hemi=ifelse(df$centlat>0,"NH","SH")
df$Sector=sector_df(df)
#converted longitude degrees (from 0->360 to 180W->180E)
df$minlon_c=lon_convert(df$minlon)
df$maxlon_c=lon_convert(df$maxlon)
df$centlon_c=lon_convert(df$centlon)



#m<-map("world",fill=TRUE,col="grey")
m<-map("world",mar=rep(.5,4))
k <- 9
my.cols <- rev(brewer.pal(k, "YlGn"))
z <- kde2d(df$centlon_c, df$centlat, n=50)
#points(df$centlon_c,df$centlat, xlab="lon", ylab="lat", pch=19, cex=.4)
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE,lwd=2)
rect(140,25,180,75,border="red")
rect(-180,25,-100,75,border="red")
#rect(40,25,140,75,border="blue")
rect(-80,25,40,75,border="red")
rect(-60,-75,30,-25,border="red")
rect(130,-75,180,-25,border="red")
rect(-180,-75,-60,-25,border="red")

#Plot the average center latitude density
#i.e. per latitude, what is center latitude density?


sec=c("ATL","PAC")
hemi=c("NH","SH")
ylim=c(0,2500)
par(mfrow=c(2,2))
for (h in hemi){
  if (h=="NH"){
    breaks=seq(24,82,2)
  }else{
    breaks=seq(-82,-24,2)

  }
  xlim=range(breaks)
  for (s in sec){
    mt=paste("Frequency of DJF center latitude values in",h,s)
    hist(df[df$Sector==s&df$Hemi==h,"centlat"],
         main=mt,xlab="centlat",
         xlim=xlim,ylim=ylim,breaks=breaks)
  }
}
par(mfrow=c(1,1))




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




# sector_blob<-function(lon0,hemi){
#   if (hemi == "NH"){
#     if ((lon0 < 40.) | (lon0 > 280.)){
#       sec = "ATL"
#     }else if ((lon0 > 140.) & (lon0 < 260.)){
#       sec = "PAC"
#     }else if ((lon0>=260.) & (lon0<=280.)){
#       sec = "CONT1"
#     }else{
#       sec = "CONT2"
#     }
#   }
#   else if (hemi == "SH"){
#     if ((lon0 < 30.) | (lon0 > 300.)){
#       sec = "ATL"
#     }else if ((lon0 > 130.) & (lon0 < 300.)){
#       sec = "PAC"
#     }else{
#       sec = "IO"
#     }
#   }
#   return(sec)
# }

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