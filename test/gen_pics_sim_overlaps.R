#load("~/block_r_data/stats_cold_spell.RData")


#library(knitr)
library(maps)
library(maptools)
library(Rmisc)
library(spatstat)
library(ggplot2)


vec_pvz<-c()
vec_pvg<-c()
vec_zg<-c()
vec_all<-c()
time_overlaps<-c()

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


brks_z<-seq(4500,6100,50)
brks_pvanom<-seq(-5.5e-6,5.5e-6,0.5e-6)
labs_pvanom<-c("-5.5",rep("",10),"0",rep("",10),"5.5")

hgt.cols<-colorRampPalette(c("purple","blue","cyan4","green","yellow", "orange","red","darkred"))(length(brks_z))
pv.cols<-colorRampPalette(c("darkred","red","white","blue","navyblue"))(length(brks_pvanom))
for (region in c("NA","NC","NP","SA","SI","SP")){
#for (region in c("NA")){
  load(sprintf('~/block_r_data/new_data/ERA_data/ERA_%s_pv_z_ghg_block_data.RData',region))
  #time_format1<-time_format

  time_hr1<-time_format
  time_hrs<-rep(seq(0,18,6),length(time_format)/4)
  load(sprintf('~/block_r_data/new_data/ERA_data/ERA_%s_z_inst_data.RData',region))
  time_hr2<-time_format
  #Need to reverse the latitude axis for each of the data
  lats_seq<-lats_seq[length(lats_seq):1]
  z_inst<-z_inst[,length(lats_seq):1,]
  z_anom<-z_anom[,length(lats_seq):1,]
  pv_anom<-pv_anom[,length(lats_seq):1,]
  ghg<-ghg[,length(lats_seq):1,]
  
  if (region=="NA" | region=="SA"){
    lons_plot=ifelse(lons_seq>179,lons_seq-360,lons_seq)
    mapr="world"
  }else{
    lons_plot=lons_seq
    mapr="world2"
  }

  for (season in c("MAM","JJA","SON","DJF")){
  #for (season in c("MAM")){
       # load(sprintf('~/block_r_data/stats_stitch_%s_%s_table.RData',season,region))

	     load(sprintf('~/block_r_data/new_data/ERA_data/ERA_stats_merged_%s_%s_table.RData',season,region))
  
    #df_tot<-df_tot_nostitch
    #df_ref<-df_tot_nostitch
    #df_tot$date_hr<-sprintf("%s_%02d",df_tot$date,df_tot$hr)
    #df_ref$date_hr<-sprintf("%s_%02d",df_ref$date,df_ref$hr)
    #print(sprintf("Length of date hr for total frame is %d, for ref is %d",nrow(df_tot),nrow(df_ref)))
    v_type<-data.frame(PV=numeric(),Z=numeric(),GHG=numeric())
    v_type_date<-data.frame(PV=numeric(),Z=numeric(),GHG=numeric(),date=character())
    for (t in sort(unique(df_tot$datehr))[1:10]){

      t1=which(time_hr1==t)
      t2=which(time_hr2==t)
      print(t1)
      print(t2)
      time_overlaps<-c(time_overlaps,sprintf("%s %s %02dZ",region, time_format1[t1],time_hrs[t1]))

      df_sub<-df_tot[df_tot$datehr==t,c("var","date","hr",
                                         "minlat","maxlat",
                                         "minlon","maxlon",
                                         "minlon_c","maxlon_c")]
      if (region=="NA" | region == "SA"){
        minlon_check<- "minlon_c"
        maxlon_check<-"maxlon_c"
      }else{
        minlon_check<- "minlon"
        maxlon_check<-"maxlon"
      }
      
      #Calculate the similarities for all three
      print("beginning to calculate similarities")
      sim_pvz<-similarity_contours(pv_anom[,,t1],z_anom[,,t1])
      vec_pvz<-c(vec_pvz,sim_pvz)
      sim_pvg<-similarity_contours(pv_anom[,,t1],ghg[,,t1])
      vec_pvg<-c(vec_pvg,sim_pvg)
      sim_zg<-similarity_contours(z_anom[,,t1],ghg[,,t1])
      vec_zg<-c(vec_zg,sim_zg)
      sim_all<-similarity_contours_3(pv_anom[,,t1],z_anom[,,t1],ghg[,,t1])
      vec_all<-c(vec_all,sim_all)
      print(c(sim_pvz,sim_pvg,sim_zg,sim_all))
      print("finished calculating similarities")
      namevec<-unique(df_sub$var)
      #Case 1: only one type of block
      if (length(namevec)==1){
        v_dir=sprintf("%s_ONLY",namevec[1])
				dir.create(file.path("~/block_r_data/figs/", v_dir), showWarnings = FALSE)
        fname<-sprintf("~/block_r_data/figs/%s/%s_%s_%02dZ_noT.png",v_dir,region,time_format1[t1],time_hrs[t1])
        print(fname)
        # png(fname,height=600,width=800)
        # print(sprintf("making picture for %d",t1))
        # map(mapr,xlim=range(lons_plot),ylim=range(lats_seq),fill=TRUE,col="grey")
        # 
        # title(sprintf("%s %02dZ, PV* (green) Z* (blue) ZG (purple)\n PV*&Z*: %f Z*&ZG: %f PV*&ZG: %f all: %f",
        #               time_format1[t1],time_hrs[t1], sim_pvz, sim_zg, sim_pvg, sim_all))
        # 
        # map.axes()
        # contour(lons_plot,lats_seq,z_inst[,,t2],levels=seq(4500,6100,50),drawlabels=TRUE,add=TRUE,col=hgt.cols,lwd=2)
        # contour(lons_plot,lats_seq,pv_anom[,,t1],levels=c(0,1),add=TRUE,col="chartreuse4",drawlabels=FALSE,lwd=5)
        # contour(lons_plot,lats_seq,z_anom[,,t1],levels=c(0,1),add=TRUE,col="cornflowerblue",drawlabels=FALSE,lwd=5)
        # contour(lons_plot,lats_seq,ghg[,,t1],levels=c(0,1),add=TRUE,col="purple",drawlabels=FALSE,lwd=5)
        # print("dev off")
        # dev.off()
      }
      #Case 2: more than one type
      else{
        #Is there overlap?
        #Case 2.1: Overlap PV&Z
        v_dir<-c()
        #Check each one of the ones in the namevec
        for (n in namevec){
          if (n=="PV" & (sim_pvg + sim_pvz<0.000001)){
            v_dir<-c(v_dir,"PV_ONLY")
          }
          if (n=="Z" & (sim_zg + sim_pvz<0.000001)){
            v_dir<-c(v_dir,"Z_ONLY")
          }
          if(n=="ZG" & (sim_pvg + sim_zg<0.000001)){
            v_dir<-c(v_dir,"ZG_ONLY")
          }
        }
        if (sim_pvz>0){
          v_dir<-c(v_dir,"PV_Z")
        }
        if(sim_pvg>0){
          v_dir<-c(v_dir,"PV_ZG")
        }
        if (sim_zg>0){
          v_dir<-c(v_dir,"Z_ZG")
        }
        if (sim_all>0){
          v_dir<-c(v_dir,"PV_Z_ZG")
        }
				dir.create(file.path("~/block_r_data/figs/", v_dir), showWarnings = FALSE)
        fname<-sprintf("~/block_r_data/figs/%s/%s_%s_%02dZ_noT.png",v_dir,region,time_format1[t1],time_hrs[t1])
        print(fname)
        # png(fname,height=600,width=800)
        # print(sprintf("making picture for %d",t1))
        # map(mapr,xlim=range(lons_plot),ylim=range(lats_seq),fill=TRUE,col="grey")
        # 
        # title(sprintf("%s %02dZ, PV* (green) Z* (blue) ZG (purple)\n PV*&Z*: %f Z*&ZG: %f PV*&ZG: %f all: %f",
        #               time_format1[t1],time_hrs[t1], sim_pvz, sim_zg, sim_pvg, sim_all))
        # 
        # map.axes()
        # contour(lons_plot,lats_seq,z_inst[,,t2],levels=seq(4500,6100,50),drawlabels=TRUE,add=TRUE,col=hgt.cols,lwd=2)
        # contour(lons_plot,lats_seq,pv_anom[,,t1],levels=c(0,1),add=TRUE,col="chartreuse4",drawlabels=FALSE,lwd=5)
        # contour(lons_plot,lats_seq,z_anom[,,t1],levels=c(0,1),add=TRUE,col="cornflowerblue",drawlabels=FALSE,lwd=5)
        # contour(lons_plot,lats_seq,ghg[,,t1],levels=c(0,1),add=TRUE,col="purple",drawlabels=FALSE,lwd=5)
        # print("dev off")
        # dev.off()

      }

      
    }


  }

}
