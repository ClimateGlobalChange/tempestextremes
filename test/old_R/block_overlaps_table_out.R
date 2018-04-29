#load("~/block_r_data/stats_cold_spell.RData")


library(knitr)
library(maps)
library(maptools)
library(Rmisc)
library(spatstat)

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

hgt.cols<-colorRampPalette(c("purple","blue","cyan4","green","yellow", "orange","red","darkred"))(34)
for (region in c("NA","NC","NP","SA","SI","SP")){
#for (region in c("NA")){
  load(sprintf('~/block_r_data/%s_pv_z_ghg_block_data.RData',region))
  time_format1<-time_format
  #time_hrs<-rep(seq(0,18,6),length(time_format1)/4)
  time_hr1<-sprintf("%s_%02d",time_format1,time_hrs)
  load(sprintf('~/block_r_data/%s_pv_z_inst_data.RData',region))
  time_hr2<-sprintf("%s_%02d",time_format2,time_hrs2)
  #Need to reverse the latitude axis for each of the data
  lats_seq<-lats_seq[length(lats_seq):1]
  z_inst<-z_inst[,length(lats_seq):1,]
  z_anom<-z_anom[,length(lats_seq):1,]
  pv_anom<-pv_anom[,length(lats_seq):1,]
  ghg<-ghg[,length(lats_seq):1,]
  
  if (region=="NA" | region=="SA"){
    lons_plot=lons_seq_c
    mapr="world"
  }else{
    lons_plot=lons_seq
    mapr="world2"
  }

  for (season in c("MAM","JJA","SON","DJF")){
  #for (season in c("MAM")){

    load(sprintf('~/block_r_data/stats_merged_%s_%s_table.RData',season,region))
    df_tot$date_hr<-sprintf("%s_%02d",df_tot$date,df_tot$hr)
    
    v_type<-data.frame(PV=numeric(),Z=numeric(),GHG=numeric())
    v_type_date<-data.frame(PV=numeric(),Z=numeric(),GHG=numeric(),date=character())
    for (t in sort(unique(df_tot$date_hr))){
      print(sprintf("Doing calculations for time %s",t))
      t1=which(time_hr1==t)
      t2=which(time_hr2==t)
      #print(t1)
      #print(t2)
      #print(time_hr1[t1])
      #print(time_hr2[t2])
      time_overlaps<-c(time_overlaps,sprintf("%s %s %02dZ",region, time_format1[t1],time_hrs[t1]))
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
      df_sub<-df_tot[df_tot$date_hr==t,c("var","date","hr",
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
      namevec<-unique(df_sub$var)
      if (length(namevec)==1){
        print(namevec)
        #print(sprintf("There is one type of block: %s",namevec[1]))
        if (namevec[1]=="GHG"){
          v_dir<-"ZG_ONLY"
        }else{
          v_dir<-sprintf("%s_ONLY",namevec[1])
        }

        fname<-sprintf("~/pics_test_case2/%s/%s_%s_%02dZ_noT.png",v_dir,region,time_format1[t1],time_hrs[t1])
        print(fname)
        
        png(fname,height=600,width=800)
        print(sprintf("making picture for %d",t1))
        map(mapr,xlim=range(lons_plot),ylim=range(lats_seq),fill=TRUE,col="grey")
        
        title(sprintf("%s %02dZ, PV* (green) Z* (blue) ZG (purple)\n PV*&Z*: %f Z*&ZG: %f PV*&ZG: %f all: %f",
                      time_format1[t1],time_hrs[t1], sim_pvz, sim_zg, sim_pvg, sim_all))
        
        map.axes()
        contour(lons_plot,lats_seq,z_inst[,,t2],levels=seq(4500,6100,50),drawlabels=TRUE,add=TRUE,col=hgt.cols,lwd=2)
        contour(lons_plot,lats_seq,pv_anom[,,t1],levels=c(0,1),add=TRUE,col="chartreuse4",drawlabels=FALSE,lwd=5)
        contour(lons_plot,lats_seq,z_anom[,,t1],levels=c(0,1),add=TRUE,col="cornflowerblue",drawlabels=FALSE,lwd=5)
        contour(lons_plot,lats_seq,ghg[,,t1],levels=c(0,1),add=TRUE,col="purple",drawlabels=FALSE,lwd=5)
        print("dev off")
        dev.off()
      }
      if (length(namevec)>1){
        #print("More than one type of block")
       # print(namevec)
        v1<-df_sub[df_sub$var==namevec[1],]
        l1<-dim(v1)[1]
        varA<-as.character(namevec[1])
        v2<-df_sub[df_sub$var==namevec[2],]
        l2<-dim(v2)[1]
        varB<-as.character(namevec[2])
        #print(sprintf("Variable A is %s and B is %s",varA,varB))
        if (length(namevec)==2){
          v3<-NULL
          varC<-NULL
    
        }
        else{
          v3<-df_sub[df_sub$var==namevec[3],]
          l3<-dim(v3)[1]
          varC<-as.character(namevec[3])
          #print(sprintf("Variable C is %s",varC))
        }
        #holds the pairs
        AB<-data.frame(A=numeric(),B=numeric())
        AC<-NULL
        BC<-NULL
        ABC<-NULL
        #New frame to hold copy of A,B,C
        AllBlobs<-data.frame(type=character(),num=numeric(),stringsAsFactors = FALSE)
        AllBlobs<-rbind(
          data.frame(type=rep("A",l1),num=seq(1,l1)),
          data.frame(type=rep("B",l2),num=seq(1,l2))
        )
        X<-data.frame(type=character(),num=numeric(), stringsAsFactors = FALSE)
    
        if (!is.null(v3)){
          AC<-data.frame(A=numeric(),C=numeric())
          BC<-data.frame(B=numeric(),C=numeric())
          ABC<-data.frame(A=numeric(),B=numeric(),C=numeric())
          
          AllBlobs<-rbind(
            AllBlobs,
            data.frame(type=rep("C",l3),num=seq(1,l3))
          )
        }
    
        #check AB
        for (a in 1:l1){
          for (b in 1:l2){
            check_AB<-check_overlaps(v1[a,"maxlat"],v1[a,"minlat"],v1[a,minlon_check],v1[a,maxlon_check],
                                     v2[b,"maxlat"],v2[b,"minlat"],v2[b,minlon_check],v2[b,maxlon_check])
            if(check_AB==TRUE){
              #Add pair to AB
              AB=rbind(AB,data.frame(A=a,B=b))
            }
          }
        }
        #check AC and BC
        if (!is.null(v3)){
          #Check AC
          for (a in 1:l1){
            for (c in 1:l3){
              check_AC<-check_overlaps(v1[a,"maxlat"],v1[a,"minlat"],v1[a,minlon_check],v1[a,maxlon_check],
                                       v3[c,"maxlat"],v3[c,"minlat"],v3[c,minlon_check],v3[c,maxlon_check])
              if (check_AC==TRUE){
                AC=rbind(AC,data.frame(A=a,C=c))
              }
            }
          }
          #check BC
          for (b in 1:l2){
            for (c in 1:l3){
              check_BC<-check_overlaps(v3[c,"maxlat"],v3[c,"minlat"],v3[c,minlon_check],v3[c,maxlon_check],
                                       v2[b,"maxlat"],v2[b,"minlat"],v2[b,minlon_check],v2[b,maxlon_check])
              if (check_BC==TRUE){
                BC=rbind(BC,data.frame(B=b,C=c))
              }
            }
          }
          
        }

        while(dim(AllBlobs)[1]>0){
          #Generate random integer to pull out and check
          randint<-as.integer(runif(1,1,dim(AllBlobs)[1]))
          blob_check1<-AllBlobs[randint,]
    
          X<-rbind(X,blob_check1)
          #Pop out the blob that will be checked
          AllBlobs<-AllBlobs[-randint,]
          blob1<-as.character(blob_check1$type[1])
          b1num<-blob_check1$num[1]
          #Check against other two types
          if (blob1=="A"){
            #print("blob is A")
            check2<-AB
            blob2<-"B"
            check3<-BC
            blob3<-"C"
            check4<-AC
          }
          else if (blob1=="B"){
            #print("blob is B")
            check2<-AB
            blob2<-"A"
            check3<-AC
            blob3<-"C"
            check4<-BC
          }
          else{
            #print("blob is C")
            check2<-AC
            blob2<-"A"
            check3<-AB
            blob3<-"B"
            check4<-BC
          }
          #We're doing 3 checks
          #First check is to find the row where the blob matches in check2
          check2_rows<-which(check2[,blob1]==b1num)
          if (length(check2_rows)>0){
            check2_vals<-check2[check2_rows,blob2]
            #Second check is where blob matches in check3
            for (z in 1:length(check2_vals)){
              #Blob of 2nd type
              row_check2<-which(AllBlobs$type==blob2 & AllBlobs$num==check2_vals[z])
              blob_check2<-AllBlobs[row_check2,]
              X<-rbind(X,blob_check2)
              AllBlobs<-AllBlobs[-row_check2,]
              b2num<-blob_check2$num[1]
              #Check against the third type
              check3_rows<-which(check3[,blob2]==b2num)
              if (length(check3_rows)>0){
                check3_vals<-check3[check3_rows,blob3]
                for (y in 1:length(check3_vals)){
                  row_check3<-which(AllBlobs$type==blob3 & AllBlobs$num==check3_vals[z])
                  blob3_pop<-AllBlobs[row_check3,]
                  X<-rbind(X,blob3_pop)
                  AllBlobs<-AllBlobs[-row_check3,]
                }
              }
            }
          }
          #Also need to run check4 to do comparison of 1 and 3
          check4_rows<-which(check4[,blob1]==b1num)
          if (length(check4_rows)>0){
            check4_vals<-check4[check4_rows,blob3]
            for (w in 1:length(check4_vals)){
              row_check4<-which(AllBlobs$type==blob3 & AllBlobs$num==check4_vals[w])
              if (length(row_check4)>0){
                blobs4_pop<-AllBlobs[row_check4,]
                X<-rbind(X,blobs4_pop)
                AllBlobs<-AllBlobs[-row_check4,]
              }
            }
          }
          type_blobs<-unique(X$type)
          if (length(type_blobs)>1){
            #print("End result:")
            #print(X)
            v_add<-data.frame(PV=0,Z=0,GHG=0)
            for (var in type_blobs){
              v_char<-as.character(var)
              varname<-ifelse(v_char=="A",varA,ifelse(v_char=="B",varB,varC))
              v_add[1,varname]=1
            }
            #print(v_add)
            v_add_date<-data.frame(PV=numeric(),Z=numeric(),GHG=numeric(),date=character(),stringsAsFactors = FALSE)
            v_add_date[1,1:3] = v_add[1,1:3]
            v_add_date[1,4] = t
            v_type<-rbind(v_type,v_add)
            v_type_date<-rbind(v_type_date,v_add_date)
            
            vnames<-c("PV","Z","ZG")
            v_dir<-NULL
            for (v in 1:3){
              if (v_add[1,v]==1){
                v_dir<-paste(v_dir,vnames[v],sep="_")
              }
            }
            v_dir<-substr(v_dir,2,10)
            

            if (length(t1)>0 & length(t2)>0){
              fname<-sprintf("~/pics_test_case2/%s/%s_%s_%02dZ_noT.png",v_dir,region,time_format1[t1],time_hrs[t1])
              print(fname)
              
              png(fname,height=600,width=800)
              print(sprintf("making picture for %d",t1))
              map(mapr,xlim=range(lons_plot),ylim=range(lats_seq),fill=TRUE,col="grey")

               title(sprintf("%s %02dZ, PV* (green) Z* (blue) ZG (purple)\n PV*&Z*: %f Z*&ZG: %f PV*&ZG: %f all: %f",
                             time_format1[t1],time_hrs[t1], sim_pvz, sim_zg, sim_pvg, sim_all))

              map.axes()
              contour(lons_plot,lats_seq,z_inst[,,t2],levels=seq(4500,6100,50),drawlabels=TRUE,add=TRUE,col=hgt.cols,lwd=2)
              contour(lons_plot,lats_seq,pv_anom[,,t1],levels=c(0,1),add=TRUE,col="chartreuse4",drawlabels=FALSE,lwd=5)
              contour(lons_plot,lats_seq,z_anom[,,t1],levels=c(0,1),add=TRUE,col="cornflowerblue",drawlabels=FALSE,lwd=5)
              contour(lons_plot,lats_seq,ghg[,,t1],levels=c(0,1),add=TRUE,col="purple",drawlabels=FALSE,lwd=5)
              print("dev off")
              dev.off()
            }
            else{
              print("No match")
            }
            
            
          }

          X<-data.frame(type=character(),num=numeric(), stringsAsFactors = FALSE)
    
        }
    
      }
      
    }


  }

}

df_vec<-data.frame(cbind(time_overlaps,vec_pvz,vec_pvg,vec_zg,vec_all))
colnames(df_vec)<-c("datetime","PV_Z","PV_ZG","Z_ZG","ALL")
hi<-0.7
lo<-0.05
file_pvz<-file("~/Desktop/sim_pvz_hi.txt")
sim_hi<-which(vec_pvz>hi)
write.table(cbind(time_overlaps[sim_hi],vec_pvz[sim_hi]),
            file=file_pvz,row.names = FALSE,col.names=FALSE,quote=FALSE)

file_pvz<-file("~/Desktop/sim_pvz_lo.txt")
sim_lo<-which(vec_pvz>0 & vec_pvz<lo)
write.table(cbind(time_overlaps[sim_lo],vec_pvz[sim_lo]),
            file=file_pvz,row.names=FALSE,col.names=FALSE,quote=FALSE)


file_pvg<-file("~/Desktop/sim_pvg_hi.txt")
sim_hi<-which(vec_pvg>hi)
write.table(cbind(time_overlaps[sim_hi],vec_pvg[sim_hi]),
            file=file_pvg,row.names=FALSE,col.names = FALSE,quote = FALSE)

file_pvg<-file("~/Desktop/sim_pvg_lo.txt")
sim_lo<-which(vec_pvg>0 & vec_pvg<lo)
write.table(cbind(time_overlaps[sim_lo],vec_pvg[sim_lo]),file=file_pvg,
            row.names = FALSE,col.names=FALSE, quote=FALSE)

file_zg<-file("~/Desktop/sim_zzg_hi.txt")
sim_hi<-which(vec_zg>hi)
write.table(cbind(time_overlaps[sim_hi],vec_zg[sim_hi]),file=file_zg,
            row.names = FALSE,col.names=FALSE, quote=FALSE)

file_zg<-file("~/Desktop/sim_zzg_lo.txt")
sim_lo<-which(vec_zg>0 & vec_zg<lo)
write.table(cbind(time_overlaps[sim_lo],vec_zg[sim_lo]),file=file_zg,
      row.names = FALSE,col.names=FALSE, quote=FALSE)

file_all<-file("~/Desktop/sim_all_hi.txt")
sim_hi<-which(vec_all>0.6)
write.table(cbind(time_overlaps[sim_hi],vec_all[sim_hi]),file=file_all,
      row.names = FALSE,col.names=FALSE, quote=FALSE)

file_all<-file("~/Desktop/sim_all_lo.txt")
sim_lo<-which(vec_all>0 & vec_all<0.05)
write.table(cbind(time_overlaps[sim_lo],vec_all[sim_lo]),file=file_all,
            row.names = FALSE,col.names=FALSE, quote=FALSE)

save(list=c("vec_pvz","vec_pvg","vec_zg","vec_all","time_overlaps"),file="~/block_r_data/similarity_vectors.RData")
