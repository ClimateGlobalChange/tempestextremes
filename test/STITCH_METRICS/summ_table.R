#This file creates a summary table using BlobStats outputs
#Input depends upon the variables that were measured

#Minlon, maxlon, minlat, maxlat, centlat, centlon, area etc

deg2rad<-function(deg) {
  return(deg * (pi/180))
}
getDistanceFromLatLonInKm<-function(lat1,lon1,lat2,lon2) {
  R = 6371; 
  dLat = deg2rad(lat2-lat1)  
  dLon = deg2rad(lon2-lon1) 
  a = sin(dLat/2) * sin(dLat/2) +
    cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * 
    sin(dLon/2) * sin(dLon/2)
  
  b = 2 * atan2(sqrt(a), sqrt(1-a))
  d = R * b 
  return(d)
}

#bcount<-1
nline<-1
gen_summary_table<-function(df_in,rfn="",textfn="",csvfn=""){
  df_summ<-data.frame(NULL)
  #print("Beginning analysis")
  #The output of these tables will be contingent upon whichever variables are available
  for (f in unique(df_in$file)){
    dsub<-df_in[df_in$file==f,]
    for (b in unique(dsub$bnum)){
      dsub2<-df_in[(df_in$file==f & df_in$bnum==b),]
      sline<-dsub2[1,]
      eline<-dsub2[nrow(dsub2),]
      df_summ[nline,"startdate"]<-sline[1,"datehour"]
      df_summ[nline,"enddate"]<-eline[1,"datehour"]
      df_summ[nline,"duration_days"]<-as.numeric(difftime(eline[1,"datehour"],
                                                          sline[1,"datehour"]))
      if (!is.null(dsub2$centlat) & !is.null(dsub2$centlon)){
        df_summ[nline,"start_centlat"]<-sline[1,"centlat"]
        df_summ[nline,"start_centlon"]<-sline[1,"centlon"]
        df_summ[nline,"end_centlat"]<-eline[1,"centlat"]
        df_summ[nline,"end_centlon"]<-eline[1,"centlon"]
        df_summ[nline,"dist_km"]<-getDistanceFromLatLonInKm(sline[1,"centlat"],
                                                            sline[1,"centlon"],
                                                            eline[1,"centlat"],
                                                            eline[1,"centlon"])
        avg_clat<-mean(dsub2$centlat)
        df_summ[nline,"zonal_dist_km"]<-getDistanceFromLatLonInKm(avg_clat,
                                                                  sline[1,"centlon"],
                                                                  avg_clat,
                                                                  eline[1,"centlon"])
        # #For zonal distance, D^2-Lat^2 = Lon ^2
        # latdist<-getDistanceFromLatLonInKm(sline[1,"centlat"],
        #                                    0,
        #                                    eline[1,"centlat"],
        #                                    0)
        # df_summ[nline,"zonal_dist_km"]<-sqrt((df_summ[nline,"dist_km"])^2-latdist^2)
        df_summ[nline,"zonal_speed_kph"]<-df_summ[nline,"zonal_dist_km"]/(df_summ[nline,"duration_days"]*24)
      }
      if (!is.null(dsub2$area_km)){
        df_summ[nline,"min_area"]<-min(dsub2$area_km)
        df_summ[nline,"max_area"]<-max(dsub2$area_km)
        df_summ[nline,"avg_area"]<-mean(dsub$area_km)
      }
      if (!is.null(dsub2$var)){
        df_summ[nline,"var"]<-dsub2[1,"var"]
      }
      df_summ[nline,"bnum"]<-b
      df_summ[nline,"file"]<-f
      #df_summ[nline,"bnum2"]<-bcount
      #bcount<-bcount+1
      nline<-nline+1
    }
  }
  if (rfn!=""){
    save(list=c("df_summ"),file=rfn) 
  }
  if (textfn!=""){
    write.table(df_summ,file=textfn,sep="\t",row.names=FALSE,quote=FALSE)
  }
  if (csvfn!=""){
    write.csv(df_summ,file=csvfn,row.names=FALSE,quote=FALSE)
  }
  #print("Returning summary table")
  return(df_summ)
}


