library("RNetCDF")
library("maps")
library("maptools")
library("mapproj")
library("abind")
library("RColorBrewer")


library("MASS")

#NOTE!!!!!!!!!!!!!!!!!!!
#NETCDF VARIABLES ARE INDEXED AS LON,LAT,TIME!!
calc_newpt<-function(x1,y1,x2,y2,xcalc){
  dy=y2-y1
  #convert points 
  nx1=lon_convert(x1)
  nx2=lon_convert(x2)
  dx=nx2-nx1
  
  slope=dy/dx
  y=slope*(xcalc-nx1) + y1
  return(y)
}

plot_pac<-function(dat,xvar,yvar,lcol="purple",lty=1,lwd=2){
  
  #check for periodic condition!
  start_i <-c(1)
  end_i<-c()
  
  for (i in 1:(nrow(dat)-1)){
    distx=abs(dat[i,xvar]-dat[i+1,xvar])
    if (distx>180){
      end_i=c(end_i,i)
      start_i=c(start_i,i+1)
    }
  }
  end_i=c(end_i,nrow(dat))
  
  for (j in 1:length(start_i)){
    
    subdat=dat[start_i[j]:end_i[j],]
    slen=nrow(subdat)
    xvec=subdat[,xvar]
    yvec=subdat[,yvar]
    #Need to connect segments
    #for all but the last segment, add row to end
    if (j<length(start_i)){
      #Look at the next point
      lastx=subdat[nrow(subdat),xvar]
      lasty=subdat[nrow(subdat),yvar]
      nextx=dat[end_i[j]+1,xvar]
      nexty=dat[end_i[j]+1,yvar]
      inty=calc_newpt(lastx,lasty,nextx,nexty,0)
      #Going eastward
      if (lastx>nextx){
        xvec=c(xvec,360)
      }else{
        xvec=c(xvec,0)
      }
      yvec=c(yvec,inty)
    }
    
    #for all but first segment, add row to beginning
    if (j>1){
      nextx=subdat[1,xvar]
      nexty=subdat[1,yvar]
      lastx=dat[start_i[j]-1,xvar]
      lasty=dat[start_i[j]-1,yvar]
      inty=calc_newpt(lastx,lasty,nextx,nexty,0)
      if (lastx<nextx){
        xvec=c(360,xvec)
      }else{
        xvec=c(0,xvec)
      }
      yvec=c(inty,yvec)
    }
    lines(xvec,yvec,col=lcol,lwd=2,lty=lty)
  }
}

produce_sectormap<-function(xlim=c(-180,180),ylim=c(-90,90)){
  m<-map('world',plot=FALSE)
  map('world',col='gray90',fill=TRUE,xlim=xlim,ylim=ylim)
  map.axes()
  
  #Boxes defining sectors
  #NH ATL
  rect(lon_convert(260),25,40,75,border="blue")
  text(((lon_convert(260)+40)/2),50,"NA")
  #NH PAC
  rect(140,25,185,75,border="blue")
  text((140+180)/2,50,"NP")
  rect(-185,25,lon_convert(260),75,border="blue")
  text((-100-180)/2,50,"NP")
  #CONTINENTS
  rect(40,25,140,75,border="red")
  text(90,50,"NC")
  
  #SH ATL
  rect(-60,-75,30,-25,border="blue")
  text(-15,-50,"SA")
  
  #SH PAC
  rect(130,-75,185,-25,border="blue")
  text((130+180)/2,-50,"SP")
  rect(-185,-75,-60,-25,border="blue")
  text(-120,-50,"SP")
  
  #Indian Ocean
  rect(30,-75,130,-25,border="blue")
  text(80,-50,"SI")
  
}
produce_sectormap2<-function(xlim=c(0,360),ylim=c(-90,90)){
  m<-map('world2',plot=FALSE)
  map('world2',xlim=xlim,ylim=ylim,fill=TRUE,col='gray90')
  map.axes()
  
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
      n1=get_index(lon,259.)
      n2=get_index(lon,40.) 
    }else if (sec=="PAC"){
      n1=get_index(lon,140.)
      n2=get_index(lon,260.)
    }else if (sec=="CONT"){
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


#-180 to 180 -> 
lon_convert2<-function(lon){
  return(ifelse(lon<0,360+lon,lon))
}


sector_df<-function(df){
  ifelse(
    df$Hemi=="NH",
    ifelse(
      #ATLANTIC
      ((df$centlon < 40.) | (df$centlon > 260.)),
      "ATL",
      #PACIFIC
      ifelse(
        ((df$centlon > 140.) & (df$centlon <= 260.)),
        "PAC",
          "CONT"
      )
    ),#SH SECTORS
    ifelse(
      #ATLANTIC
      ((df$centlon < 30.) | (df$centlon > 300.)),
      "ATL",
      ifelse(
        #PACIFIC
        ((df$centlon >= 130.) & (df$centlon <= 300.)),
        "PAC",
        #INDIAN OCEAN
        "IO"
      )
    )
  )
}