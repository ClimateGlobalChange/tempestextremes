library("RNetCDF")
library("maps")
library("maptools")
library("mapproj")
library("abind")
library("RColorBrewer")


library("MASS")

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