args=commandArgs(trailingOnly=TRUE)
data=args[1]
sector=args[2]
season=args[3]
#for (data in c("ERA","MERRA")){
#for (season in c("MAM","JJA","SON","DJF")){
#  for (sector in c("NA","NC","NP","SA","SI","SP")){
#for (season in c("JJA")){
#  for (sector in c("NP")){
    fname<-sprintf("~/block_r_data/%s_stats_merged_%s_%s_table.RData",data,season,sector)
    #load(fname)
    f2<-sprintf("~/block_r_data/%s_%s_%s_stats_stitch_table.RData",data,season,sector)
    f3<-sprintf("~/block_r_data/%s_%s_%s_stats_nostitch_table.RData",data,season,sector)
    load(f2)
    load(f3)
    df_tot_stitch$area_km<-(df_tot_stitch$area)*(4*pi*(6.371e6)^2)/1e6
    df_tot_nostitch$area_km<-(df_tot_nostitch$area)*(4*pi*(6.371e6)^2)/1e6
    df_comm<-merge(df_tot_stitch,df_tot_nostitch,by=c("date","area_km","hr","var","tstep"))
    df_comm$datehr<-sprintf("%s_%02d",df_comm$date,df_comm$hr)
    df_all<-merge(df_tot_stitch,df_tot_nostitch,by=c("date","area_km","hr","var","tstep"), all=TRUE)
    
    df_istot<-df_all[is.na(df_all$bnum.y),]
    df_istot$datehr<-sprintf("%s_%02d",df_istot$date,df_istot$hr)
    df_isnot<-df_all[is.na(df_all$bnum.x),]
    df_isnot$datehr<-sprintf("%s_%02d",df_isnot$date,df_isnot$hr)
    
    df_tot<-df_comm[,c(1:18,length(names(df_comm)))]
    df_tot$bnum2<-df_tot$bnum.x
    min1<-ifelse((sector=="NA"|sector=="SA"),"minlon_c.x","minlon.x")
    max1<-ifelse((sector=="NA"|sector=="SA"),"maxlon_c.x","maxlon.x")
    for (t in unique(df_istot$datehr)){
      df_check<-df_istot[df_istot$datehr==t,]
      df_sub<-df_isnot[df_isnot$datehr==t,c(1:5,19:length(names(df_comm)))]
      colnames(df_sub)<-gsub(".y",".x",names(df_sub))
      df_sub$bnum2<-df_sub$bnum.x
      for (n in 1:nrow(df_check)){
        clatmin<-df_check[n,"minlat.x"]
        clatmax<-df_check[n,"maxlat.x"]
        clonmin<-df_check[n,min1]
        clonmax<-df_check[n,max1]
        for (y in 1:nrow(df_sub)){
          #Does it fall within the bounds of the big stitched blob?
          blatmin<-df_sub[y,"minlat.x"]
          blatmax<-df_sub[y,"maxlat.x"]
          blonmin<-df_sub[y,min1]
          blonmax<-df_sub[y,max1]
          if (((blatmin>=clatmin)&(blatmax<=clatmax)&
               (blonmin>=clonmin)&(blonmax<=clonmax)&
               (df_check[n,"var"]==df_sub[y,"var"]))){
            #print("block falls inside original merged block")
            df_sub[y,"bnum.x"]<-df_check[n,"bnum.x"]
            df_tot<-rbind(df_tot,df_sub[y,])
            
          }else{
            print("outside of range")
          }
        }
      }
    }
    save(list=c("df_tot","df_tot_stitch","df_tot_nostitch"),file=fname)
#  }
#}
#}
