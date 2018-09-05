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
gen_summary_table<-function(df_in,keep_merge=TRUE,
                            rfn="",textfn="",csvfn="",
                            df_summ_name=""){
  df_summ<-data.frame(NULL)
  df_summ_name<-ifelse(df_summ_name=="","df_summ",df_summ_name)
  #print("Beginning analysis")
  #The output of these tables will be contingent upon whichever variables are available
  for (f in unique(df_in$file)){
    dsub<-df_in[df_in$file==f,]
    for (v in unique(dsub$var)){
      for (b in unique(dsub$bnum)){
        dsub2<-df_in[(df_in$file==f & df_in$bnum==b & df_in$var==v),]
        dsub2<-dsub2[order(dsub2$datehour),]
        merged_blob<-FALSE
        if (!is.null(dsub2$bnum2)){
          #Check if there are any instances in which bnum!=bnum2
          for (a in 1:nrow(dsub2)){
            if (dsub2[a,"bnum"]!=dsub2[a,"bnum2"]){
              merged_blob<-TRUE
              break
            }
          }
        }
        
        sline<-dsub2[1,]
        eline<-dsub2[nrow(dsub2),]
        df_summ[nline,"startdate"]<-sline[1,"datehour"]
        df_summ[nline,"enddate"]<-eline[1,"datehour"]
        diff_days<-as.numeric(difftime(as.Date(eline[1,"datehour"]),as.Date(sline[1,"datehour"]),units="days"))
        hs<-as.numeric(strftime(sline[1,"datehour"],format="%H"))
        he<-as.numeric(strftime(eline[1,"datehour"],format="%H"))
        hdiff<-(he-hs)/24
        df_summ[nline,"duration_days"]<-diff_days+hdiff
        df_summ[nline,"merged"]<-ifelse(merged_blob==FALSE,"NO","YES")
        if (!is.null(dsub2$centlat) & !is.null(dsub2$centlon)){
          df_summ[nline,"start_centlat"]<-sline[1,"centlat"]
          df_summ[nline,"start_centlon"]<-sline[1,"centlon"]
          df_summ[nline,"end_centlat"]<-eline[1,"centlat"]
          df_summ[nline,"end_centlon"]<-eline[1,"centlon"]
          df_summ[nline,"dist_km"]<-getDistanceFromLatLonInKm(sline[1,"centlat"],
                                                              sline[1,"centlon"],
                                                              eline[1,"centlat"],
                                                              eline[1,"centlon"])
          avg_clat<-(sline[1,"centlat"]+eline[1,"centlat"])*0.5
          df_summ[nline,"zonal_dist_km"]<-getDistanceFromLatLonInKm(avg_clat,
                                                                    sline[1,"centlon"],
                                                                    avg_clat,
                                                                    eline[1,"centlon"])
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

  }
  if (keep_merge==FALSE){
    df_summ<-df_summ[df_summ$merged=="NO",]
  }
  if (rfn!=""){
    assign(df_summ_name,df_summ)
    save(list=c(df_summ_name),file=rfn) 
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


