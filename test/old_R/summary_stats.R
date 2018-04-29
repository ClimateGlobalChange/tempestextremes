library(ggplot2)
library(akima)
deg2rad<-function(deg) {
  return(deg * (pi/180))
}
getDistanceFromLatLonInKm<-function(lat1,lon1,lat2,lon2) {
  R = 6371; # Radius of the earth in km
  dLat = deg2rad(lat2-lat1)  
  dLon = deg2rad(lon2-lon1) 
  a = sin(dLat/2) * sin(dLat/2) +
    cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * 
    sin(dLon/2) * sin(dLon/2)
  
  b = 2 * atan2(sqrt(a), sqrt(1-a))
  d = R * b # Distance in km
  # if (is.na(d)){
  #   print(sprintf("getting NAs for %f, %f, %f, %f",lat1,lon1,lat2,lon2))
  # }
  return(d)
}

for (season in c("MAM","JJA","SON","DJF")){
#if sector is Atlantic, use centlon_c
#else use centlon
df_summ<-data.frame(var=character(),bnum=numeric(),
                    startdate=character(),starthr=numeric(),
                    enddate=character(),endhr=numeric(),duration=numeric(),
                    minlat=numeric(),maxlat=numeric(),
                    minlon=numeric(),maxlon=numeric(),
                    minlatdiff=numeric(),maxlatdiff=numeric(),
                    minlondiff=numeric(),maxlondiff=numeric(),
                    avglatdiff=numeric(),avglondiff=numeric(),
                    startclat=numeric(),startclon=numeric(),
                    endclat=numeric(),endclon=numeric(),
                    minclat=numeric(),minclon=numeric(),
                    maxclat=numeric(),maxclon=numeric(),
                    avgclat=numeric(),avgclon=numeric(),
                    minarea=numeric(),maxarea=numeric(),
                    avgarea=numeric(),dist=numeric(),sum_dist=numeric(),
                    sector=character(),avgkmh=numeric(),avgzonalmps=numeric(),
                    stringsAsFactors = FALSE)

nr=1
for (sector in c("NA","NC","NP","SA","SI","SP")){
  load(sprintf("~/block_r_data/stats_merged_%s_%s_table.RData",season,sector))
  #match the bnums from the stitched table with the total table
  test<-merge (df_tot,df_tot_stitch,by=c("date","hr","var","area"))
  for (b in sort(unique(test$bnum))){
    test_sub<-test[test$bnum==b,]
    test_sub$latdiff<-test_sub$maxlat.x-test_sub$minlat.x
    nt<-nrow(test_sub)
    varsub<-unique(as.character(test_sub$var))
    if (length(varsub)>1){
      print("Error! 2 blocks in subset")
    }else{
      t<-as.numeric(test_sub[nt,"date"]-test_sub[1,"date"])
      if(sector=="NA" | sector=="SA"){
        minlon_sec<-test_sub$minlon_c.x
        maxlon_sec<-test_sub$maxlon_c.x
        centlon_sec<-test_sub$centlon_c.x
        test_sub$londiff<-test_sub$maxlon_c.x-test_sub$minlon_c.x
      }else{
        minlon_sec<-test_sub$minlon.x
        maxlon_sec<-test_sub$maxlon.x
        centlon_sec<-test_sub$centlon.x
        test_sub$londiff<-test_sub$maxlon.x-test_sub$minlon.x
      }
      for (x in 2:nt){
        # print(sprintf("for x=%d, %f, %f, %f, %f",x,test_sub[x-1,"centlat.x"],
        #               centlon_sec[x-1],
        #               test_sub[x,"centlat.x"],
        #               centlon_sec[x]))
        test_sub[x,"dist_traveled"]<-getDistanceFromLatLonInKm(test_sub[x-1,"centlat.x"],
                                                               centlon_sec[x-1],
                                                               test_sub[x,"centlat.x"],
                                                               centlon_sec[x])
        #Zonal dist in km
        avglattest<-(test_sub[x-1,"centlat.x"]+test_sub[x,"centlat.x"])*0.5
        test_sub[x,"zonal_dist"]<-getDistanceFromLatLonInKm(avglattest,
                                                               centlon_sec[x-1],
                                                               avglattest,
                                                               centlon_sec[x])
        test_sub[x,"kmh"]<-test_sub[x,"zonal_dist"]/6
        test_sub[x,"mps"]<-test_sub[x,"zonal_dist"]*(10/(6*36))
      }
      test_sub[1,"zonal_dist"]<-0
      test_sub[1,"dist_traveled"]<-0
      test_sub[1,"mps"]<-0
      test_sub[1,"kmh"]<-0
  
      if (t>=5){
  
        if (any(is.na(centlon_sec))){
          print("getting NAs!")
          print(centlon_sec)
        }        
  
        df_summ[nr,"duration"]<-t
          df_summ[nr,"bnum"]<-b
          df_summ[nr,"var"]<-varsub
          df_summ[nr,"startdate"]<-as.character(test_sub[1,"date"])
          df_summ[nr,"enddate"]<-as.character(test_sub[nt,"date"])
          df_summ[nr,"starthr"]<-test_sub[1,"hr"]
          df_summ[nr,"endhr"]<-test_sub[nt,"hr"]
          df_summ[nr,"startclat"]<-test_sub[1,"centlat.x"]
          df_summ[nr,"endclat"]<-test_sub[nt,"centlat.x"]
          df_summ[nr,"startclon"]<-centlon_sec[1]
          df_summ[nr,"endclon"]<-centlon_sec[nt]
          df_summ[nr,"minlat"]<-min(test_sub$minlat.x)
          df_summ[nr,"maxlat"]<-max(test_sub$maxlat.x)
          df_summ[nr,"minlon"]<-min(minlon_sec)
          df_summ[nr,"maxlon"]<-max(maxlon_sec)
          df_summ[nr,"minlatdiff"]<-min(test_sub$latdiff)
          df_summ[nr,"maxlatdiff"]<-max(test_sub$latdiff)
          df_summ[nr,"minlondiff"]<-min(test_sub$londiff)
          df_summ[nr,"maxlondiff"]<-max(test_sub$londiff)
          df_summ[nr,"avglatdiff"]<-mean(test_sub$latdiff)
          df_summ[nr,"avglondiff"]<-mean(test_sub$londiff)
          df_summ[nr,"minclat"]<-min(test_sub$centlat.x)
          df_summ[nr,"minclon"]<-min(centlon_sec)
          df_summ[nr,"maxclat"]<-max(test_sub$centlat.x)
          df_summ[nr,"maxclon"]<-max(centlon_sec)
          df_summ[nr,"avgclat"]<-mean(test_sub$centlat.x)
          df_summ[nr,"avgclon"]<-mean(centlon_sec)
          df_summ[nr,"minarea"]<-min(test_sub$area)
          df_summ[nr,"maxarea"]<-max(test_sub$area)
          df_summ[nr,"avgarea"]<-mean(test_sub$area)
          df_summ[nr,"dist"]<-getDistanceFromLatLonInKm(test_sub[1,"centlat.x"],
                                                        centlon_sec[1],
                                                        test_sub[nt,"centlat.x"],
                                                        centlon_sec[nt])
          df_summ[nr,"sum_dist"]<-sum(test_sub$dist_traveled)
          df_summ[nr,"sector"]<-sector
          df_summ[nr,"avgzonalmps"]<-mean(test_sub$mps)
          df_summ[nr,"avgkmh"]<-mean(test_sub$kmh)
          # difflon<-abs(df_summ[nr,"maxlon"]-df_summ[nr,"minlon"])
          # difflon_c<-abs(df_summ[nr,"maxlon_c"]-df_summ[nr,"minlon_c"])
          # if (difflon<difflon_c){
          #   print("using difflon")
          # }else{
          #   print("using difflon_c")
          # }
          nr<-nr+1
      }
    }
  }
}
df_summ$avgclon_c<-ifelse(df_summ$avgclon<0,df_summ$avgclon+360,df_summ$avgclon)
save(list=c("df_summ"),file=sprintf("~/block_r_data/%s_summ_stats.RData",season))
}

breaks_duration<-seq(5,40,5)
duration_cols<-colorRampPalette(c("green","blue",
                                  "purple","red"))(length(breaks_duration))
gridfld<-interp(df_summ$avgclon,df_summ$avgclat,df_summ$duration)
df_grid<-as.data.frame(interp2xyz(gridfld))
# df_summ$avgclat<-round(df_summ$avgclat,digits=2)
# df_summ$avgclon<-round(df_summ$avgclon,digits=2)



mn<-ifelse((sector=="NA"|sector=="SA"),"world","world2")
m<-map_data(mn)

# ggplot(data=test)+
#   coord_cartesian(xlim=c(0,360),
#                   ylim=c(-90,90),expand=FALSE) +
#   geom_map(data=m,map=m,aes(map_id=region))+
#   geom_path(data=subset(test,area>1200000),
#             aes(x=centlon.x,y=centlat.x,color=var,group=bnum))


ggplot(data=df_summ,expand=FALSE)+
  facet_wrap(~sector)+
  geom_freqpoly(aes(x=dist,y=..count..,color=var),binwidth=1000)

ggplot(data=df_summ,expand=FALSE)+
  facet_wrap(~sector)+
  geom_freqpoly(aes(x=avgarea,y=..count..,color=var),binwidth=100000)


##### WHICH METHOD?

ggplot(data=df_summ,expand=FALSE)+
  facet_wrap(~sector)+
  geom_freqpoly(aes(x=duration,y=..count..,color=var),binwidth=5)+
  xlim(0,50)

ggplot(data=df_summ,expand=FALSE)+
  facet_wrap(~sector)+
  geom_density(aes(x=duration,color=var))+
  ylim(0,.4)+xlim(0,50)

ggplot(data=df_summ)+ 
  coord_cartesian(xlim=c(0,360),
                  ylim=c(-90,90),expand=FALSE) +
  geom_map(data= m, map = m, aes(map_id=region))+
  scale_size_continuous(breaks=seq(0,80,20))+
  geom_point(data=subset(df_summ,duration>20),
             aes(x=avgclon_c,y=avgclat,color=var,size=duration),alpha=0.75)


ggplot(data=df_summ)+ 
  coord_cartesian(xlim=c(0,360),
                  ylim=c(-90,90),expand=FALSE) +
  geom_map(data= m, map = m, aes(map_id=region))+
  geom_point(data=subset(df_summ,duration<6),
             aes(x=avgclon_c,y=avgclat,color=var),alpha=0.75)

######
gplot<-ggplot(data=df_summ)+ 
  coord_cartesian(xlim=c(0,360),
                  ylim=c(-90,90),expand=FALSE) +
  geom_map(data= m, map = m, aes(map_id=region))+
  geom_point(aes(x=avgclon_c,y=avgclat,color=var,size=dist),alpha=0.5)+
  scale_size_continuous(breaks=seq(5,35,5))
  # geom_point(data=subset(df_summ,duration<6),
  #             aes(x=avgclon_c,y=avgclat,color=var,size=duration),alpha=0.75)
  # scale_color_gradientn(colors=duration_cols) +
  # geom_contour(data=df_grid,aes(x=x,y=y,z=z))