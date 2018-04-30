library(abind)
library(reshape2)
library(ggplot2)
library(spatstat)
check_overlaps<-function(Alt,Alb,All,Alr,
                         Blt,Blb,Bll,Blr){
  #Note: lons must deal with periodic boundary
  #use centered lons for Atlantic
  #Check that the latitudes overlap
  if ((Alt<Blb) | (Alb>Blt)){
    ###print("lats don't overlap")
    return(FALSE)
  }
  #Check that the longitudes overlap
  else if((Alr<Bll) | (All>Blr)){
    ###print("lons don't overlap")
    return(FALSE)
  }else{
    ###print("overlap!")
    return(TRUE)
  }
}

similarity_contours<-function(arr1,arr2){
  arr1[which(arr1>0)]<-1
  arr2[which(arr2>0)]<-1
  sum1<-sum(arr1)
  sum2<-sum(arr2)
  if (sum1>0 & sum2>0){
    arr_sum<-arr1+arr2
    arr_sum[which(arr_sum>0)]<-1
    arr.im<-as.im(arr_sum)
    arr.connect<-connected(arr.im,background=0)
    n_lev<-length(levels(arr.connect))
    #What are the relative sizes of the blocks?
    cluster_points<-c()
    sim_vec<-c()
    if (n_lev>0){
      for (l in levels(arr.connect$v)){
        arr1_sub<-arr1[which(arr.connect$v==l)]
        arr2_sub<-arr2[which(arr.connect$v==l)]
        #How many common 1's are there?
        arr_connect_sum<-arr1_sub + arr2_sub
        length_both<-length(arr_connect_sum)
        n_both<-length(which(arr_connect_sum>1))
        #save size of cluster
        cluster_points<-c(cluster_points,length_both)
        #find similarity
        frac_both<-n_both/length_both
        sim_vec<-c(sim_vec,frac_both)
      }
      #Find the fraction of points 
      ntotalpts<-sum(cluster_points)
      cluster_points<-cluster_points/ntotalpts
      #weight the similarities
      weighted_sim<-cluster_points*sim_vec
      #weighted
      sim_total<-sum(weighted_sim)
      #unweighted
      sim_unweight<-mean(sim_vec)
    }
  }else{
    sim_total<-0
    sim_unweight<-0
  }
  return(sim_total)
}
similarity_contours_3<-function(arr1,arr2,arr3){
  arr1[which(arr1>0)]<-1
  arr2[which(arr2>0)]<-1
  arr3[which(arr3>0)]<-1
  sum1<-sum(arr1)
  sum2<-sum(arr2)
  sum3<-sum(arr3)
  if (sum1>0 & sum2>0 & sum3>0){
    arr_sum<-arr1+arr2+arr3
    arr_sum[which(arr_sum>0)]<-1
    arr.im<-as.im(arr_sum)
    arr.connect<-connected(arr.im,background=0)
    n_lev<-length(levels(arr.connect))
    #What are the relative sizes of the blocks?
    cluster_points<-c()
    sim_vec<-c()
    if (n_lev>0){
      for (l in levels(arr.connect$v)){
        arr1_sub<-arr1[which(arr.connect$v==l)]
        arr2_sub<-arr2[which(arr.connect$v==l)]
        arr3_sub<-arr3[which(arr.connect$v==l)]
        #How many common 1's are there?
        arr_connect_sum<-arr1_sub + arr2_sub + arr3_sub
        length_all<-length(arr_connect_sum)
        n_all<-length(which(arr_connect_sum>2))
        #save size of cluster
        cluster_points<-c(cluster_points,length_all)
        #find similarity
        frac_all<-n_all/length_all
        sim_vec<-c(sim_vec,frac_all)
      }
      #Find the fraction of points 
      ntotalpts<-sum(cluster_points)
      cluster_points<-cluster_points/ntotalpts
      #weight the similarities
      weighted_sim<-cluster_points*sim_vec
      #weighted
      sim_total<-sum(weighted_sim)
      #unweighted
      sim_unweight<-mean(sim_vec)
    }
  }else{
    sim_total<-0
    sim_unweight<-0
  }
  return(sim_total)
}

prob_vec<-function(df,dfo,lat="ALL"){
  
  if (lat=="ALL"){
    dfp<-df
    dfz<-df
    dfg<-df
    dfo<-dfo
  }else{
    dfp<-df[df$var=="PV" & df$HILO==lat,]
    dfz<-df[df$var=="Z" & df$HILO==lat,]
    dfg<-df[df$var=="GHG" & df$HILO==lat,]
    dfo<-dfo[(!is.na(dfo$PVZHILO==lat) | 
                !is.na(dfo$PVGHGHILO==lat) | 
                !is.na(dfo$ZGHGHILO==lat)),]
  }
  PVcount<-c(nrow(dfp[dfp$var=="PV",]))
  Zcount<-c(nrow(dfz[dfz$var=="Z",]))
  GHGcount<-c(nrow(dfg[dfg$var=="GHG",]))
  PVZcount<-0
  PVGcount<-0
  ZZGcount<-0
  ALLcount<-0
  
  for (d in sort(unique(dfo$datehr))){
    #How many unique overlaps are there per time step?
    dsub<-dfo[dfo$datehr==d,]
    if (lat=="ALL"){
      nPV<-length(unique(dsub[!is.na(dsub$PVbnum2),"PVbnum2"]))
      nZ<-length(unique(dsub[!is.na(dsub$Zbnum2),"Zbnum2"]))
      nGHG<-length(unique(dsub[!is.na(dsub$GHGbnum2),"GHGbnum2"])) 
      nPVZ<-length(unique(dsub[!is.na(dsub$PV_Z),"PV_Z"]))
      nZZG<-length(unique(dsub[!is.na(dsub$Z_GHG),"Z_GHG"]))
      nPVZG<-length(unique(dsub[!is.na(dsub$PV_GHG),"PV_GHG"]))
      nALL<-length(unique(dsub[!is.na(dsub$ALL),"ALL"]))
    }
    else{
      nPV<-length(unique(dsub[!is.na(dsub$PVbnum2) & dsub$PVHILO==lat,"PVbnum2"]))
      nZ<-length(unique(dsub[!is.na(dsub$Zbnum2) & dsub$ZHILO==lat,"Zbnum2"]))
      nGHG<-length(unique(dsub[!is.na(dsub$GHGbnum2) & dsub$GHGHILO==lat,"GHGbnum2"])) 
      nPVZ<-length(unique(dsub[!is.na(dsub$PV_Z) & dsub$PVZHILO==lat,"PV_Z"]))
      nZZG<-length(unique(dsub[!is.na(dsub$Z_GHG) & dsub$ZGHGHILO==lat,"Z_GHG"]))
      nPVZG<-length(unique(dsub[!is.na(dsub$PV_GHG) & dsub$PVGHGHILO==lat,"PV_GHG"]))
      nALL<-length(unique(dsub[!is.na(dsub$ALL),"ALL"]))
    }

    
    PVZ<-NULL
    PVZG<-NULL
    ZZG<-NULL
    ALL<-NULL
    if (nPV>0 & nZ>0 & nPVZ>0){
      PVZ<-min(nPVZ,nPV,nZ)
      PVZcount<-PVZcount+PVZ
    }
    if (nPV>0 & nGHG>0 & nPVZG>0){
      PVZG<-min(nPV,nGHG,nPVZG)
      PVGcount<-PVGcount+PVZG
    }
    if (nZ>0 & nGHG>0 & nZZG>0){
      ZZG<-min(nZ, nGHG, nZZG)
      ZZGcount<-ZZGcount+ZZG
    }
    if (nZ>0 & nGHG>0 & nPV>0 & nALL>0){
      ALL<-min(nZ, nGHG, nPV,nALL)
      ALLcount<-ALLcount+ALL
    }
    
  }
  
  #Probabilities
  PVZgivenZ<-PVZcount/Zcount
  PVGgivenG<-PVGcount/GHGcount
  PVZgivenP<-PVZcount/PVcount
  ZZGgivenG<-ZZGcount/GHGcount
  PVGgivenP<-PVGcount/PVcount
  ZZGgivenZ<-ZZGcount/Zcount
  
  pvec<-c(PVZgivenZ,PVZgivenP,PVGgivenG,PVGgivenP,ZZGgivenG,ZZGgivenZ)
  return(pvec)
}


for (region in c("NA","NC","NP","SA","SI","SP")){
#for (region in c("SA","SI","SP")){
  load(sprintf("~/block_r_data/%s_pv_z_ghg_block_data.RData",region))
  
  lat_plot<-rev(lats_seq)
  
  if (region=="NA" | region=="SA"){
    lon_plot<-lons_seq_c
    clon_var<-"centlon_c.x"
    lon_max<-"maxlon_c.x"
    lon_min<-"minlon_c.x"
    lon_c<-"centlon_c.x"
  }else{
    lon_plot<-lons_seq
    clon_var<-"centlon.x"
    lon_max<-"maxlon.x"
    lon_min<-"minlon.x"
    lon_c<-"centlon.x"
  }
  for (season in c("MAM","JJA","SON","DJF")){
    load(sprintf("~/block_r_data/stats_merged_%s_%s_table.RData",season,region))
    df_tot$HILO<-ifelse(abs(df_tot$centlat.x)<40,"LOLAT","HILAT")
    # load(sprintf("~/block_r_data/%s_summ_stats.RData",season))
    # 
    # df_summ$KMPH<-df_summ$dist/(df_summ$duration*24)
    # df_summ$HILO<-ifelse(abs(df_summ$avgclat)<40,"LOWLAT",ifelse(abs(df_summ$avgclat)>60,"HILAT","MIDLAT"))
    # df_summ$HILO<-factor(df_summ$HILO,levels=c("HILAT","MIDLAT","LOWLAT"))
    
    df_overlaps<-data.frame(datehr=character(),PVbnum=numeric(),Zbnum=numeric(),GHGbnum=numeric(),
                            PVbnum2=numeric(),Zbnum2=numeric(),GHGbnum2=numeric(),
                            PVminlat=numeric(),PVminlon=numeric(),
                            PVmaxlat=numeric(),PVmaxlon=numeric(),
                            PVcentlat=numeric(),PVcentlon=numeric(),
                            Zminlat=numeric(),Zminlon=numeric(),
                            Zmaxlat=numeric(),Zmaxlon=numeric(),
                            Zcentlat=numeric(),Zcentlon=numeric(),
                            GHGminlat=numeric(),GHGminlon=numeric(),
                            GHGmaxlat=numeric(),GHGmaxlon=numeric(),
                            GHGcentlat=numeric(),GHGcentlon=numeric(),
                            PV_Z=numeric(),PV_GHG=numeric(),Z_GHG=numeric(),ALL=numeric(),
                            PV_Zsize=numeric(),PV_GHGsize=numeric(),Z_GHGsize=numeric(),ALLsize=numeric(),
                            stringsAsFactors = FALSE
    )
    nr=1
    time_vec<-sprintf("%s_%02d",time_format,time_hrs)
    for (d in sort(unique(df_tot$datehr))){
      #print(sprintf("Date %s",d))
      ncount<-1
      overcopy<-df_overlaps[1,]
      overcopy[,]<-NA
      #Data subset for time step
      dsub=df_tot[df_tot$datehr==d,]
      dsub<-dsub[!duplicated(dsub),]
      di<-which(time_vec==d)
      #If only one type of block, add and skip calcs
      nvar<-length(unique(dsub$var))
      varname<-as.character(sort(unique(as.character(dsub$var))))
      if (nvar==1){
        #print("Only one type of block")
        # for (n in 1:nrow(dsub)){
        #   df_overlaps[nr,"datehr"]<-d
        #   df_overlaps[nr,sprintf("%sminlat",varname)]<-dsub[n,"minlat.x"]
        #   df_overlaps[nr,sprintf("%smaxlat",varname)]<-dsub[n,"maxlat.x"]
        #   df_overlaps[nr,sprintf("%sminlon",varname)]<-dsub[n,lon_min]
        #   df_overlaps[nr,sprintf("%smaxlon",varname)]<-dsub[n,lon_max]
        #   df_overlaps[nr,sprintf("%sbnum",varname)]<-dsub[n,"bnum.x"]
        #   df_overlaps[nr,sprintf("%sbnum2",varname)]<-dsub[n,"bnum2"]
        # 
        #   nr<-nr+1
        #}
      }else{
        #A bit more complicated! Check overlaps each against each
        
        if (varname[1]=="GHG"){
          v1<-ghg
          df1<-dsub[dsub$var=="GHG",]
          if (varname[2]=="PV"){
            v2<-pv_anom
            df2<-dsub[dsub$var=="PV",]
            name_12<-"PV_GHG"
          }else{
            v2<-z_anom
            df2<-dsub[dsub$var=="Z",]
            name_12<-"Z_GHG"
          }
        }else{
          v1<-pv_anom
          df1<-dsub[dsub$var=="PV",]
          v2<-z_anom
          df2<-dsub[dsub$var=="Z",]
          name_12<-"PV_Z"
        }
        if (length(varname)>2){
          v3<-z_anom
          v3val<-100
          df3<-dsub[dsub$var=="Z",]
          name_23<-"PV_Z"
          name_13<-"Z_GHG"
        }else{
          v3<-NULL
          df3<-NULL
          name_23<-NULL
          name_13<-NULL
          v3count<-NULL
        }
        
        #Loop through the variables
        for (b1 in 1:nrow(df1)){
          it1<-which(lats_seq==as.integer(df1[b1,"maxlat.x"]))
          ib1<-which(lats_seq==as.integer(df1[b1,"minlat.x"]))
          il1<-which(lon_plot==as.integer(df1[b1,lon_min]))
          ir1<-which(lon_plot==as.integer(df1[b1,lon_max]))
          v1slice<-v1[,,di]
          v1slice[c(1:min(il1,ir1),max(il1,ir1):dim(v1)[1]),]<-0
          v1slice[,c(1:min(it1,ib1),max(it1,ib1):dim(v1)[2])]<-0
          v1count<-length(which(v1slice>0))
          for (b2 in 1:nrow(df2)){
            it2<-which(lats_seq==as.integer(df2[b2,"maxlat.x"]))
            ib2<-which(lats_seq==as.integer(df2[b2,"minlat.x"]))
            il2<-which(lon_plot==as.integer(df2[b2,lon_min]))
            ir2<-which(lon_plot==as.integer(df2[b2,lon_max]))
            v2slice<-v2[,,di]
            v2slice[c(1:min(il2,ir2),max(il2,ir2):dim(v2)[1]),]<-0
            v2slice[,c(1:min(it2,ib2),max(it2,ib2):dim(v2)[2])]<-0
            v2count<-length(which(v2slice>0))
            #Check overlap between 1 and 2
            s12<-NULL
            s13<-NULL
            s23<-NULL
            sAll<-NULL
            v3slice<-NULL
            v13count<-NULL
            v23count<-NULL
            vAllcount<-NULL
            
            s12<-similarity_contours(v1slice,v2slice)
            v12count<-length(which((v1slice + v2slice)>0))
            if (!is.null(v3)){
              #print("There are 3 variables")
              for (b3 in 1:nrow(df3)){
                #print("CHeck 3: inside null")
                it3<-which(lats_seq==as.integer(df3[b3,"maxlat.x"]))
                ib3<-which(lats_seq==as.integer(df3[b3,"minlat.x"]))
                il3<-which(lon_plot==as.integer(df3[b3,lon_min]))
                ir3<-which(lon_plot==as.integer(df3[b3,lon_max]))
                v3slice<-v3[,,di]
                v3slice[c(1:min(il3,ir3),max(il3,ir3):dim(v3)[1]),]<-0
                v3slice[,c(1:min(it3,ib3),max(it3,ib3):dim(v3)[2])]<-0
                v3count<-length(which(v3slice>0))
                
                v13count<-length(which((v1slice + v3slice)>0))
                v23count<-length(which((v3slice + v2slice)>0))
                vAllcount<-length(which((v1slice+v2slice+v3slice)>0))
                overcopy[ncount,name_12]<-s12
                #Check overlap between 1 and 3
                s13<-similarity_contours(v1slice,v3slice)
                #Check overlap between 2 and 3
                s23<-similarity_contours(v2slice,v3slice)
                #Check overlap between all 3
                sAll<-similarity_contours_3(v1slice,v2slice,v3slice)
                overcopy[ncount,"datehr"]<-d
                overcopy[ncount,sprintf("%sminlat",varname[1])]<-df1[b1,"minlat.x"]
                overcopy[ncount,sprintf("%smaxlat",varname[1])]<-df1[b1,"maxlat.x"]
                overcopy[ncount,sprintf("%sminlon",varname[1])]<-df1[b1,lon_min]
                overcopy[ncount,sprintf("%smaxlon",varname[1])]<-df1[b1,lon_max]
                overcopy[ncount,sprintf("%sbnum",varname[1])]<-df1[b1,"bnum.x"]
                overcopy[ncount,sprintf("%sbnum2",varname[1])]<-df1[b1,"bnum2"]
                overcopy[ncount,sprintf("%scentlat",varname[1])]<-df1[b1,"centlat.x"]
                overcopy[ncount,sprintf("%scentlon",varname[1])]<-df1[b1,lon_c]
                overcopy[ncount,sprintf("%sminlat",varname[2])]<-df2[b2,"minlat.x"]
                overcopy[ncount,sprintf("%smaxlat",varname[2])]<-df2[b2,"maxlat.x"]
                overcopy[ncount,sprintf("%sminlon",varname[2])]<-df2[b2,lon_min]
                overcopy[ncount,sprintf("%smaxlon",varname[2])]<-df2[b2,lon_max]
                overcopy[ncount,sprintf("%sbnum",varname[2])]<-df2[b2,"bnum.x"]
                overcopy[ncount,sprintf("%sbnum2",varname[2])]<-df2[b2,"bnum2"]
                overcopy[ncount,sprintf("%scentlat",varname[2])]<-df2[b2,"centlat.x"]
                overcopy[ncount,sprintf("%scentlon",varname[2])]<-df2[b2,lon_c]
                overcopy[ncount,sprintf("%sminlat",varname[3])]<-df3[b3,"minlat.x"]
                overcopy[ncount,sprintf("%smaxlat",varname[3])]<-df3[b3,"maxlat.x"]
                overcopy[ncount,sprintf("%sminlon",varname[3])]<-df3[b3,lon_min]
                overcopy[ncount,sprintf("%smaxlon",varname[3])]<-df3[b3,lon_max]
                overcopy[ncount,sprintf("%sbnum",varname[3])]<-df3[b3,"bnum.x"]
                overcopy[ncount,sprintf("%sbnum2",varname[3])]<-df3[b3,"bnum2"]
                overcopy[ncount,sprintf("%scentlat",varname[3])]<-df3[b3,"centlat.x"]
                overcopy[ncount,sprintf("%scentlon",varname[3])]<-df3[b3,lon_c]
                overcopy[ncount,"ALL"]<-sAll
                overcopy[ncount,name_13]<-s13
                overcopy[ncount,name_23]<-s23
                overcopy[ncount,sprintf("%ssize",varname[1])]<-v1count
                overcopy[ncount,sprintf("%ssize",varname[2])]<-v2count
                overcopy[ncount,sprintf("%ssize",varname[3])]<-v3count
                overcopy[ncount,sprintf("%ssize",name_12)]<-v12count
                overcopy[ncount,sprintf("%ssize",name_13)]<-v13count
                overcopy[ncount,sprintf("%ssize",name_23)]<-v23count
                overcopy[ncount,"ALLsize"]<-vAllcount
                #print("At the end of adding to overcopy")
                ncount<-ncount+1
              }
              
            }else{
              #print("There are 2 variables")
              overcopy[ncount,name_12]<-s12
              #Check overlap between 1 and 3
              s13<-similarity_contours(v1slice,v3slice)
              #Check overlap between 2 and 3
              s23<-similarity_contours(v2slice,v3slice)
              #Check overlap between all 3
              sAll<-similarity_contours_3(v1slice,v2slice,v3slice)
              overcopy[ncount,"datehr"]<-d
              overcopy[ncount,sprintf("%sminlat",varname[1])]<-df1[b1,"minlat.x"]
              overcopy[ncount,sprintf("%smaxlat",varname[1])]<-df1[b1,"maxlat.x"]
              overcopy[ncount,sprintf("%sminlon",varname[1])]<-df1[b1,lon_min]
              overcopy[ncount,sprintf("%smaxlon",varname[1])]<-df1[b1,lon_max]
              overcopy[ncount,sprintf("%sbnum",varname[1])]<-df1[b1,"bnum.x"]
              overcopy[ncount,sprintf("%sbnum2",varname[1])]<-df1[b1,"bnum2"]
              overcopy[ncount,sprintf("%scentlat",varname[1])]<-df1[b1,"centlat.x"]
              overcopy[ncount,sprintf("%scentlon",varname[1])]<-df1[b1,lon_c]
              #print("filled v1")
              overcopy[ncount,sprintf("%sminlat",varname[2])]<-df2[b2,"minlat.x"]
              overcopy[ncount,sprintf("%smaxlat",varname[2])]<-df2[b2,"maxlat.x"]
              overcopy[ncount,sprintf("%sminlon",varname[2])]<-df2[b2,lon_min]
              overcopy[ncount,sprintf("%smaxlon",varname[2])]<-df2[b2,lon_max]
              overcopy[ncount,sprintf("%sbnum",varname[2])]<-df2[b2,"bnum.x"]
              overcopy[ncount,sprintf("%sbnum2",varname[2])]<-df2[b2,"bnum2"]
              overcopy[ncount,sprintf("%scentlat",varname[2])]<-df2[b2,"centlat.x"]
              overcopy[ncount,sprintf("%scentlon",varname[2])]<-df2[b2,lon_c]
              #print("filled v2")
              overcopy[ncount,sprintf("%ssize",varname[1])]<-v1count
              overcopy[ncount,sprintf("%ssize",varname[2])]<-v2count
              overcopy[ncount,sprintf("%ssize",name_12)]<-v12count
              
             # print("At the end of adding to overcopy")
              ncount<-ncount+1
            }
            
          }
          
        }
        #Unique combos of vars
        #print("Beginning melt")
        df_meltpz<-melt(overcopy[,c("PVbnum2","Zbnum2","PV_Z")],
                        id.vars<-c("PVbnum2","Zbnum2"))
        df_meltpz_nonzero<-df_meltpz[(df_meltpz$value>0 & 
                                        !is.na(df_meltpz$value)),c(1:2,4)] 
        df_pz<-df_meltpz_nonzero[!duplicated(df_meltpz_nonzero),]
        df_meltzg<-melt(overcopy[,c("GHGbnum2","Zbnum2","Z_GHG")],
                        id.vars<-c("GHGbnum2","Zbnum2"))
        df_meltzg_nonzero<-df_meltzg[(df_meltzg$value>0 & 
                                        !is.na(df_meltzg$value)),c(1:2,4)]
        df_zg<-df_meltzg_nonzero[!duplicated(df_meltzg_nonzero),]
        df_meltpg<-melt(overcopy[,c("PVbnum2","GHGbnum2","PV_GHG")],
                        id.vars<-c("PVbnum2","GHGbnum2"))
        df_meltpg_nonzero<-df_meltpg[(df_meltpg$value>0 & 
                                        !is.na(df_meltpg$value)),c(1:2,4)]
        df_pg<-df_meltpg_nonzero[!duplicated(df_meltpg_nonzero),]
        df_meltall<-melt(overcopy[,c("PVbnum2","Zbnum2","GHGbnum2","ALL")],
                         id.vars<-c("PVbnum2","Zbnum2","GHGbnum2"))
        df_meltall_nonzero<-df_meltall[(df_meltall$value>0 & 
                                          !is.na(df_meltall$value)),c(1:3,5)]
        df_all<-df_meltall_nonzero[!duplicated(df_meltall_nonzero),]
        
        df_merge1<-merge(df_all[,1:3],df_pz[,1:2],id.vars=c("PVbnum2","Zbnum2"),all=T)
        df_merge2<-merge(df_merge1,df_pg[,1:2],id.vars=c("PVbnum2","GHGbnum2"),all=T)
        df_merge3<-merge(df_merge2,df_zg[,1:2],id.vars=c("GHGbnum2","Zbnum2"),all=T)
        
        if (nrow(df_merge3)>0){
          for (r in 1:nrow(df_merge3)){
            df_overlaps[nr,"datehr"]<-d
            bp<-df_merge3[r,"PVbnum2"]
            bg<-df_merge3[r,"GHGbnum2"]
            bz<-df_merge3[r,"Zbnum2"]
            
            for (v in varname){
              #Match the bnum to the proper row
              col<-which(colnames(df_merge3)==sprintf("%sbnum2",v))
              #Blob number
              bv<-df_merge3[r,col]
              #Match blob number to blob size
              # brow<-overcopy[overcopy[,sprintf("%sbnum",v)]==bv,]
              # bsize<-unique(brow[,sprintf("%ssize",v)])
              if (!is.na(bv)){
                ri<-which(dsub$var==v & dsub$bnum2==bv)
                df_overlaps[nr,sprintf("%sminlat",v)]<-dsub[ri,"minlat.x"]
                df_overlaps[nr,sprintf("%smaxlat",v)]<-dsub[ri,"maxlat.x"]
                df_overlaps[nr,sprintf("%sminlon",v)]<-dsub[ri,lon_min]
                df_overlaps[nr,sprintf("%smaxlon",v)]<-dsub[ri,lon_max]
                df_overlaps[nr,sprintf("%sbnum",v)]<-dsub[ri,"bnum.x"]
                df_overlaps[nr,sprintf("%sbnum2",v)]<-dsub[ri,"bnum2"]
                df_overlaps[nr,sprintf("%scentlat",v)]<-dsub[ri,"centlat.x"]
                df_overlaps[nr,sprintf("%scentlon",v)]<-dsub[ri,lon_c]
              }
              if (!is.na(bp)){
                sizep<-unique(overcopy[overcopy$PVbnum2==bp,"PVsize"])
                df_overlaps[nr,"PVsize"]<-sizep
                if (!is.na(bg)){
                  bpg<-unique(overcopy[overcopy$PVbnum2==bp & 
                                         overcopy$GHGbnum2==bg,"PV_GHG"])
                  sizepg<-unique(overcopy[overcopy$PVbnum2==bp & 
                                            overcopy$GHGbnum2==bg,"PV_GHGsize"])
                  sizeg<-unique(overcopy[overcopy$GHGbnum2==bg,"GHGsize"])
                  df_overlaps[nr,"GHGsize"]<-sizeg   
                  if (length(bpg)>0){
                    bpg<-bpg[!is.na(bpg)]
                  }
                  df_overlaps[nr,"PV_GHG"]<-bpg
                  df_overlaps[nr,"PV_GHGsize"]<-sizepg      
                }
                if (!is.na(bz)){
                  bpz<-unique(overcopy[overcopy$PVbnum2==bp & 
                                         overcopy$Zbnum2==bz,"PV_Z"])
                  sizepz<-unique(overcopy[overcopy$PVbnum2==bp & 
                                            overcopy$Zbnum2==bz,"PV_Zsize"])
                  sizez<-unique(overcopy[overcopy$Zbnum2==bz,"Zsize"])
                  df_overlaps[nr,"Zsize"]<-sizez  
                  if (length(bpz)>0){
                    bpz<-bpz[!is.na(bpz)]
                  }
                  df_overlaps[nr,"PV_Z"]<-bpz
                  df_overlaps[nr,"PV_Zsize"]<-sizepz
                }
              }
              if (!is.na(bg) & !is.na(bz)){
                bzg<-unique(overcopy[overcopy$GHGbnum2==bg & 
                                       overcopy$Zbnum2==bz,"Z_GHG"])
                
                sizezg<-unique(overcopy[overcopy$GHGbnum2==bg & 
                                          overcopy$Zbnum2==bz,"Z_GHGsize"])
                
                sizeg<-unique(overcopy[overcopy$GHGbnum2==bg,"GHGsize"])
                sizez<-unique(overcopy[overcopy$Zbnum2==bz,"Zsize"])
                
                df_overlaps[nr,"GHGsize"]<-sizeg
                df_overlaps[nr,"Zsize"]<-sizez  
                if (length(bzg)>0){
                  bzg<-bzg[!is.na(bzg)]
                }
                df_overlaps[nr,"Z_GHG"]<-bzg
                df_overlaps[nr,"Z_GHGsize"]<-sizezg
                if (!is.na(bp)){
                  ball<-unique(overcopy[overcopy$PVbnum2==bp & 
                                          overcopy$Zbnum2==bz &
                                          overcopy$GHGbnum2==bg,"ALL"])
                  sizeall<-unique(overcopy[overcopy$PVbnum2==bp & 
                                             overcopy$Zbnum2==bz &
                                             overcopy$GHGbnum2==bg,"ALLsize"])
                  if (length(ball)>0){
                    ball<-ball[!is.na(ball)]
                  }
                  df_overlaps[nr,"ALL"]<-ball
                  df_overlaps[nr,"ALLsize"]<-sizeall
                }
              }
            }
            nr<-nr+1
          }
        }else{
          for (r in 1:nrow(dsub)){
            v<-dsub[r,"var"]
            bv<-dsub[r,sprintf("bnum2",v)]
            brow<-overcopy[overcopy[,sprintf("%sbnum2",v)]==bv,]
            bsize<-unique(brow[,sprintf("%ssize",v)])
            df_overlaps[nr,"datehr"]<-d
            df_overlaps[nr,sprintf("%sminlat",v)]<-dsub[r,"minlat.x"]
            df_overlaps[nr,sprintf("%smaxlat",v)]<-dsub[r,"maxlat.x"]
            df_overlaps[nr,sprintf("%sminlon",v)]<-dsub[r,lon_min]
            df_overlaps[nr,sprintf("%smaxlon",v)]<-dsub[r,lon_max]
            df_overlaps[nr,sprintf("%sbnum",v)]<-dsub[r,"bnum.x"]
            df_overlaps[nr,sprintf("%sbnum2",v)]<-dsub[r,"bnum2"]
            df_overlaps[nr,sprintf("%scentlat",v)]<-dsub[r,"centlat.x"]
            df_overlaps[nr,sprintf("%scentlon",v)]<-dsub[r,lon_c]
            df_overlaps[nr,sprintf("%ssize",v)]<-bsize
            nr<-nr+1
          }
        }
      }
    }
    
    df_overlaps$PVHILO<-ifelse(is.na(df_overlaps$PVcentlat),NA,
                                   ifelse(abs(df_overlaps$PVcentlat)<40,"LOLAT","HILAT"))
    df_overlaps$ZHILO<-ifelse(is.na(df_overlaps$Zcentlat),NA,
                                  ifelse(abs(df_overlaps$Zcentlat)<40,"LOLAT","HILAT"))
    df_overlaps$GHGHILO<-ifelse(is.na(df_overlaps$GHGcentlat),NA,
                                    ifelse(abs(df_overlaps$GHGcentlat)<40,"LOLAT","HILAT"))
    df_overlaps$PVZlat<-(df_overlaps$PVcentlat + df_overlaps$Zcentlat)*0.5
    df_overlaps$PVGlat<-(df_overlaps$PVcentlat + df_overlaps$GHGcentlat)*0.5
    df_overlaps$ZZGlat<-(df_overlaps$GHGcentlat + df_overlaps$Zcentlat)*0.5
    df_overlaps$PVZHILO<-ifelse(is.na(df_overlaps$PVZlat),NA,
                                    ifelse(abs(df_overlaps$PVZlat)<40,"LOLAT","HILAT"))
    df_overlaps$PVGHGHILO<-ifelse(is.na(df_overlaps$PVGlat),NA,
                                      ifelse(abs(df_overlaps$PVGlat)<40,"LOLAT","HILAT"))
    df_overlaps$ZGHGHILO<-ifelse(is.na(df_overlaps$ZZGlat),NA,
                                     ifelse(abs(df_overlaps$ZZGlat)<40,"LOLAT","HILAT"))
    
    
    # We have: Similarity
    # Relative block sizes
    # For overlap, unique per day
    pvec_all<-prob_vec(df_tot,df_overlaps)
    pvec_hi<-prob_vec(df_tot,df_overlaps,lat="HILAT")    
    pvec_lo<-prob_vec(df_tot,df_overlaps,lat="LOLAT")
        #print(sprintf("Probs for %s %s: %f", region, season, pvec))
    
    
    
    #Save this info 
    save(list=c("df_tot","df_overlaps","pvec_all","pvec_hi","pvec_lo"),
         file=sprintf("~/block_r_data/%s_%s_prob_overlaps_sim.RData",season,region))
  }
}


#load(sprintf("~/block_r_data/%s_pv_z_inst_anom_data.RData",region))
#load(sprintf("~/block_r_data/%s_pv_z_inst_data.RData",region))

#df_summ$SUMMKMPH<-df_summ$sum_dist/(df_summ$duration*24)


#Start times!
#Is there some sort of lag between start time of one vs another? 
#i.e. do they 

# mmin<-ifelse(season=="MAM",3,ifelse(season=="JJA",6,9))
# mmax<-ifelse(season=="MAM",5,ifelse(season=="JJA",8,11))

# latmin<-ifelse(substr(region,1,1)=="N",25,-75)
# latmax<-ifelse(substr(region,1,1)=="N",75,-25)



#New dataframe with similarity comparisons and overlaps
#Goals for today: GET SOME WRITING DONE



#Find the start date
#Find the end date






# ## Richard doesn't like this. Save for reference ----
# if (season=="DJF"){
#   dateinds<-which((as.numeric(format(time_format,"%m"))>11 | as.numeric(format(time_format,"%m"))<3))
# }else{
#   dateinds<-which((as.numeric(format(time_format,"%m"))>=mmin & as.numeric(format(time_format,"%m"))<=mmax))
#   
# }
# 
# pv_sub<-pv_anom[,,dateinds]
# pv_sub[which(pv_sub>0)]<-1
# pv_avg<-apply(pv_sub,c(1,2),mean)
# pv_avg[which(pv_avg<0.02)]<-NA
# z_sub<-z_anom[,,dateinds]
# z_sub[which(z_sub>1)]<-1
# z_avg<-apply(z_sub,c(1,2),mean)
# z_avg[which(z_avg<0.02)]<-NA
# zg_sub<-ghg[,,dateinds]
# zg_sub[which(zg_sub>1)]<-1
# zg_avg<-apply(zg_sub,c(1,2),mean)
# zg_avg[which(zg_avg<0.02)]<-NA
# all_avg<-abind(pv_avg,z_avg,zg_avg,along=3)
# 
# #Version 2 of similarity composites
# #Step 1: Find the max density of lat/lon coordinates. This will be the center of the box
# #For each grid point, find the mean of all three densities. 
# #Subtract the difference between the highest and the lowest
# 
# avg_dens<-apply(all_avg,c(1,2),mean)
# max_dens<-apply(all_avg,c(1,2),max)
# min_dens<-apply(all_avg,c(1,2),min)
# 
# diff_dens<-avg_dens-(max_dens-min_dens)*.33
# maxind<-which(diff_dens==max(diff_dens,na.rm = T),arr.ind = T)
# latcoord<-lats_seq[maxind[,2]]
# loncoord<-lon_plot[maxind[,1]]
# 
# contour(lon_plot,lat_plot,pv_avg[,ll:1],levels=seq(.0,.2,.01),col="darkgreen")
# title(sprintf("Avg frequencies at max point: PV=%f, Z=%f, ZG=%f",
#               pv_avg[maxind],z_avg[maxind],zg_avg[maxind]))
# contour(lon_plot,lat_plot,z_avg[,ll:1],levels=seq(.05,.2,.01),col="blue",add=T)
# contour(lon_plot,lat_plot,zg_avg[,ll:1],levels=seq(.0,.2,.01),col="purple",add=T)
# contour(lon_plot,lat_plot,diff_dens[,ll:1],col="red",levels=c(0.01,0.06,.01),add=T,lwd=3)
# 
# points(x=loncoord,y=latcoord,cex=1.5,pch=16)
# 
# latboxmax<-latcoord+4
# latboxmin<-latcoord-5
# lonboxmax<-loncoord+7
# lonboxmin<-loncoord-7
# 
# rect(lonboxmin,latboxmin,lonboxmax,latboxmax,lwd=3,border="blue")
# 
# df_tot_sub<-df_tot[df_tot[,"centlat.x"]>latboxmin &
#                      df_tot[,"centlat.x"]<latboxmax &
#                      df_tot[,clon_var]>lonboxmin &
#                      df_tot[,clon_var]<lonboxmax &
#                      df_tot[,"maxlat.x"]<latmax &
#                      df_tot[,"minlat.x"]>latmin &
#                      df_tot[,lon_max]<max(lon_plot) &
#                      df_tot[,lon_min]>min(lon_plot),]
# 
# #Similarity between fields at these times?
# 
# df_pv<-df_tot_sub[df_tot_sub$var=="PV",]
# df_z<-df_tot_sub[df_tot_sub$var=="Z",]
# df_zg<-df_tot_sub[df_tot_sub$var=="GHG",]
# 
# df_date_sim<-date.frame(datehr=character(),PV_Z=numeric(),PV_ZG=numeric(),
#                         Z_ZG=numeric(),ALL=numeric())
# 
# datehr_vec<-sprintf("%s_%s",time_format,time_hrs)
# nr<-1
# for (d in sort(unique(df_tot_sub$datehr))[1:20]){
#   df_sub_sub<-df_tot_sub[df_tot_sub$datehr==d,]
#   print(df_sub_sub)
#   #Get the date index
#   dind=which(datehr_vec==d)
#   #Get the max extent of the blocks being compared 
#   #(in order to isolate just the blocks of interest)
#   
# }
# 
# 
# 
# # For NA DJF, there are 23 PV blocks, 18 Z blocks, 25 ZG blocks
# # For the number of cases, choose the type with the least number of unique blocks
# lp<-length(unique(df_pv$bnum.x))
# lz<-length(unique(df_z$bnum.x))
# lg<-length(unique(df_zg$bnum.x))
# 
# minl<-min(c(lp,lz,lg))
# 
# varmin<-ifelse(minl==lp,"PV",ifelse(minl==lz,"Z","GHG"))
# 
# 
