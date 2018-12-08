for (sector in c("NA","NC","NP","SA","SI","SP")){
  if (sector=="NA" | sector == "SA"){
    minlon_check<-"minlon_c.y"
    minlon_test<-"minlon_c.x"
    maxlon_check<-"maxlon_c.y"
    maxlon_test<-"maxlon_c.x"
  }
  else{
    minlon_check<-"minlon.y"
    minlon_test<-"minlon.x"
    maxlon_check<-"maxlon.y"
    maxlon_test<-"maxlon.x"
  }
  for (season in c("MAM","JJA","SON","DJF")){
    load(sprintf("~/block_r_data/stats_nostitch_%s_%s_table.RData",season,sector))
    load(sprintf("~/block_r_data/stats_stitch_%s_%s_table.RData",season,sector))
    
    
    #Note! This only gets the ones where there's 100% match
    #Doesn't get the ones where 2 blobs merge into 1
    #~8000 
    df_new<-merge(df_tot_nostitch,df_tot_stitch,by=c("date","hr","var","tstep","area"))
    df_all<-merge(df_tot_nostitch,df_tot_stitch,by=c("date","hr","var","tstep","area"),all=TRUE)
    df_all$date_hr<-sprintf("%s_%02d",df_all$date,df_all$hr)
    #This is all the nostitch data omitting the non-matching stitch data
    df_omit_stitch<-df_all[which(is.na(df_all$bnum.y)),]
    #This is all the stitch data omitting the non-matching nostitch data
    df_omit_nostitch<-df_all[which(is.na(df_all$bnum.x)),]
    
    df_col<-colnames(df_new)
    
    rcount<-nrow(df_new)+1
    #Now we need to find matches between the stitch and nostitch data
    for (var in c("GHG","PV","Z")){
      sub_stitchdata<-df_omit_nostitch[df_omit_nostitch$var==var,]
      sub_nostitch<-df_omit_stitch[df_omit_stitch$var==var,]
      for (d in unique(sub_stitchdata$date_hr)){
        sub_stitch<-sub_stitchdata[sub_stitchdata$date_hr==d,
                                   c("minlat.y","maxlat.y","minlon.y","maxlon.y",
                                     "minlon_c.y","maxlon_c.y")]
        
        sub_sub<-sub_nostitch[sub_nostitch$date_hr==d,]
        #The bounding box for stitch should contain the bounding boxes for nostitch
        for (r in 1:nrow(sub_sub)){
          sub_check<-sub_sub[r,]
          #check the lats
          if (sub_check$minlat.x[1]>=sub_stitch$minlat.y[1] & sub_check$maxlat.x[1]<=sub_stitch$maxlat.y[1]){
            #Check the longitudes
            if (sub_check[1,minlon_test]>=sub_stitch[1,minlon_check] & sub_check[1,maxlon_test] <= sub_stitch[1,maxlon_check]){
              for (col in df_col){
                df_new[rcount,col]<-sub_check[1,col]
              }
              rcount<-rcount+1
            }
          }
        }
        
      }
    }
    for (col in c("minlat.x","maxlat.x","minlon.x","maxlon.x","minlon_c.x","maxlon_c.x",
                  "centlat.x","centlon.x","centlon_c.x")){
      newname<-strsplit(col,"\\.")[[1]][1]
      colnames(df_new)[which(colnames(df_new)==col)]<-newname
    }
    df_tot<-df_new[,1:15]
    fname<-sprintf("~/block_r_data/stats_merged_%s_%s_table.RData",season,sector)
    save(list=c("df_tot","df_tot_stitch","df_tot_nostitch"),file=fname)
  }
}

