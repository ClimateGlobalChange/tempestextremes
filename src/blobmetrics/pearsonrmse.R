require(reshape2)
require(akima)

#Returns in the -180->180 range
lon_convert<-function(lon){
  distFrom180=lon-180.
  return(ifelse(
    distFrom180<0,
    lon,
    -(180-distFrom180)
  ))
}
#Returns in the 0->360 range
lon_convert2<-function(lon){
  return(ifelse(lon<0,360+lon,lon))
}

# #Finds the Pearson pattern correlation between the averaged climatologies of the two datasets
# pearson_arr<-function(arr1,arr2,lat1,lat2,lon1,lon2,interp=FALSE,centered=FALSE){
#   
#   longdata1<-melt(arr1,value.name="V1")
#   longdata1$lon<-lon1[longdata1$Var1]
#   longdata1$lat<-lat1[longdata1$Var2]
#   
#   if (interp==TRUE){
#     arr2_analyze<-linint(arr2,lat2,lon2,lat1,lon1)
#   }else{
#     arr2_analyze<-arr2
#   }
#   
#   longdata2<-melt(arr2_analyze,value.name="V2")
#   longdata2$lon<-lon1[longdata2$Var1]
#   longdata2$lat<-lat1[longdata2$Var2]
#   
#   longdata<-merge(longdata1,longdata2,by=c("lon","lat"))
#   #Create a cosine latitude column
#   longdata$coslat<-cos(longdata$lat*pi/180)
#   longdata$cosV1<-longdata$V1*longdata$coslat
#   longdata$cosV2<-longdata$V2*longdata$coslat
#   
#   Vproduct<-longdata$cosV1*longdata$cosV2
#   
#   V1bar<-mean(longdata$cosV1)
#   V2bar<-mean(longdata$cosV2)
#   
#   V1diff<-longdata$cosV1-V1bar
#   V2diff<-longdata$cosV2-V2bar
#   VdiffProduct<-V1diff*V2diff
#   
#   if (centered==TRUE){
#     r<-sum(VdiffProduct)/sqrt(sum(V1diff*V1diff)*sum(V2diff*V2diff))    
#   }else{
#     r<-sum(Vproduct)/sqrt(sum(longdata$cosV1*longdata$cosV1)*sum(longdata$cosV2*longdata$cosV2))
#   }
#   
#   return(r)
# }

pearson_permute<-function(arr1,arr2,lat1,lat2,lon1,lon2,interp=FALSE,centered=FALSE){
  N_s<-1000
  
  longdata1<-melt(arr1,value.name="V1")
  longdata1$lon<-lon1[longdata1$Var1]
  longdata1$lat<-lat1[longdata1$Var2]
  
  if (interp==TRUE){
    arr2_analyze<-linint(arr2,lat2,lon2,lat1,lon1)
  }else{
    arr2_analyze<-arr2
  }
  
  longdata2<-melt(arr2_analyze,value.name="V2")
  longdata2$lon<-lon1[longdata2$Var1]
  longdata2$lat<-lat1[longdata2$Var2]
  
  longdata<-merge(longdata1,longdata2,by=c("lon","lat"))
  #Create a cosine latitude column
  longdata$coslat<-cos(longdata$lat*pi/180)
  longdata$cosV1<-longdata$V1*longdata$coslat
  longdata$cosV2<-longdata$V2*longdata$coslat
  
  v1<-longdata$cosV1
  v2<-longdata$cosV2
  #Observed global correlation:
  cor_obs<-as.vector(cor(v1,v2))
  
  N<-nrow(longdata)
  
  cvals<-c()
  cmean<-c()
  for (x in 1:10000){
    #if N is too large, take a sample
    if (N>N_s){
      ids<-sample(1:N,N_s,replace=F)
      longdata_sample<-longdata[ids,]
    }else{
      longdata_sample<-longdata
    }
    #Randomly permute v2
    v1sample<-longdata_sample$cosV1
    v2sample<-longdata_sample$cosV2
    v2permute<-sample(v2sample,length(v2sample),replace=F)
    cor_new<-cor(v1sample,v2permute)
    
    cvals[x]<-as.vector(cor_new)
    cmean[x]<-as.vector(cor(v1sample,v2sample))
  }
  pval<-mean(abs(cvals)>=abs(cor_obs))
  vals<-list(cor.value=cor_obs,p.value=pval,c.vec=cmean)
  return(vals)
}
rmse_permute<-function(arr1,arr2,lat1,lat2,lon1,lon2,interp=FALSE){
  N_s<-1000
  longdata1<-melt(arr1,value.name="V1")
  longdata1$lon<-lon1[longdata1$Var1]
  longdata1$lat<-lat1[longdata1$Var2]
  
  if (interp==TRUE){
    arr2_analyze<-linint(arr2,lat2,lon2,lat1,lon1)
  }else{
    arr2_analyze<-arr2
  }
  
  longdata2<-melt(arr2_analyze,value.name="V2")
  longdata2$lon<-lon1[longdata2$Var1]
  longdata2$lat<-lat1[longdata2$Var2]
  
  longdata<-merge(longdata1,longdata2,by=c("lon","lat"))
  #Create a cosine latitude column
  longdata$coslat<-cos(longdata$lat*pi/180)
  longdata$cosV1<-longdata$V1*longdata$coslat
  longdata$cosV2<-longdata$V2*longdata$coslat
  
  diff_data2<-(longdata$cosV1-longdata$cosV2)^2
  rmse_obs<-sqrt(sum(diff_data2)/length(diff_data2))
  rmse_perm<-c()
  rmse_sub<-c()
  N<-nrow(longdata)  
  for (x in 1:10000){
    if (N>N_s){
      ids<-sample(1:N,N_s,replace=F)
      longdata_sample<-longdata[ids,]
    }else{
      longdata_sample<-longdata
    }
    v1sample<-longdata_sample$cosV1
    v2sample<-longdata_sample$cosV2
    v2permute<-sample(v2sample,length(v2sample),replace=F)
    diff_data2<-(v1sample-v2permute)^2
    rmse_perm[x]<-sqrt(sum(diff_data2)/length(diff_data2))
    diff_data2.new<-(v1sample-v2sample)^2
    rmse_sub[x]<-sqrt(sum(diff_data2.new)/length(diff_data2.new))
  }
  #Is the RMSE on the left or right tail?
  ncount<-ifelse(rmse_obs>=median(rmse_perm),sum(rmse_perm>=rmse_obs),sum(rmse_perm<=rmse_obs))
  pval<-ncount/length(rmse_perm)
  vals<-list(rmse.value=rmse_obs,p.value=pval,rmse.vec=rmse_sub)
   #print(sprintf("AME is %f and RMSE is %f",ame,rmse))
  return(vals)
}


# #Calculates the rmse of 2 wrt 1
# rmse_calc<-function(arr1,arr2,lat1,lat2,lon1,lon2,interp=FALSE){
#   longdata1<-melt(arr1,value.name="V1")
#   longdata1$lon<-lon1[longdata1$Var1]
#   longdata1$lat<-lat1[longdata1$Var2]
#   
#   if (interp==TRUE){
#     arr2_analyze<-linint(arr2,lat2,lon2,lat1,lon1)
#   }else{
#     arr2_analyze<-arr2
#   }
#   
#   longdata2<-melt(arr2_analyze,value.name="V2")
#   longdata2$lon<-lon1[longdata2$Var1]
#   longdata2$lat<-lat1[longdata2$Var2]
#   
#   longdata<-merge(longdata1,longdata2,by=c("lon","lat"))
#   #Create a cosine latitude column
#   longdata$coslat<-cos(longdata$lat*pi/180)
#   longdata$cosV1<-longdata$V1*longdata$coslat
#   longdata$cosV2<-longdata$V2*longdata$coslat
#   
#   diff_data2<-(longdata$cosV1-longdata$cosV2)^2
#   ame<-sum(abs(longdata$cosV1-longdata$cosV2))/nrow(longdata)
#   rmse<-sqrt(sum(diff_data2)/length(diff_data2))
#   #print(sprintf("AME is %f and RMSE is %f",ame,rmse))
#   return(rmse)
# }

pearson_rmse<-function(V1,V2,blobs1,time_format1,lat1,lon1,blobs2,time_format2,lat2,lon2,
                       regrid=FALSE,time_intersect=FALSE,
                       rfn_pr="",txt_pr=""){
  if (time_intersect==TRUE){
    timeNetCDF<-intersect(time_format1,time_format2)
    #Get the time indices that are in the NetCDF intersect for V1
    time_sub1<-which(time_format1 %in% timeNetCDF)
    #Get the time indices that are in the NetCDF intersect for V2
    time_sub2<-which(time_format2 %in% timeNetCDF)
  }else{
    time_sub1<-seq(1,length(time_format1))
    time_sub2<-seq(1,length(time_format2))
  }
  b1avg<-apply(blobs1[,,time_sub1],c(1,2),mean)
  b2avg<-apply(blobs2[,,time_sub2],c(1,2),mean)
  
  #Does the longitude axis have the same range of values?
  #Create longitude axes that have both 0->360 and -180->180 axes
  lon1_180<-lon_convert(lon1)
  lon2_180<-lon_convert(lon2)
  lon1_360<-lon_convert2(lon1)
  lon2_360<-lon_convert2(lon2)
  
  #Use 180 or 360?
  if (lon1_360[1]>lon1_360[length(lon1_360)] | lon2_360[1]>lon2_360[length(lon2_360)]){
    lon1_analyze<-lon1_180
    lon2_analyze<-lon2_180
  }else{
    lon1_analyze<-round(lon1_360,3)
    lon2_analyze<-round(lon2_360,3)
  }
  
  pearson_calc<-pearson_permute(b1avg,b2avg,round(lat1,3),round(lat2,3),
                           round(lon1_analyze,3),round(lon2_analyze,3),interp=regrid)
  pearson_num<-pearson_calc$cor.value
  pearson_p<-pearson_calc$p.value
  rmse_calc<-rmse_permute(b1avg,b2avg,round(lat1,3),round(lat2,3),
                      round(lon1_analyze,3),round(lon2_analyze,3),interp=regrid)
  rmse_num<-rmse_calc$rmse.value
  rmse_p<-rmse_calc$p.value
  saved_names<-c("V1","V2","pearson_num","rmse_num","pearson_p","rmse_p")
  if (rfn_pr!=""){
    save(list=saved_names,file=rfn_pr)
    print(sprintf("Saved %s to file",rfn_pr))
  }
  if (txt_pr!=""){
    sink(txt_pr)
    cat(sprintf("Variable 1:%s \n",V1))
    cat(sprintf("Pearson correlation between average of 1 and 2: %f \n",pearson_num))
    cat(sprintf("p-value for Pearson correlation: %f\n",pearson_p))
    cat(sprintf("RMSE between average of 1 and 2: %f \n",rmse_num))
    cat(sprintf("p-value for RMSE: %f\n",rmse_p))
    sink()
  }
}
