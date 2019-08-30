require(reshape2)
pearson_arr<-function(arr1,arr2,lat1,lat2,lon1,lon2,interp=FALSE,centered=FALSE){

  longdata1<-melt(arr1,value.name="V1")
  longdata1$lon<-lon1[longdata1$Var1]
  longdata1$lat<-lat1[longdata1$Var2]
  
  
  longdata2<-melt(arr2,value.name="V2")
  longdata2$lon<-lon2[longdata2$Var1]
  longdata2$lat<-lat2[longdata2$Var2]
  if (interp==TRUE){
    
    temp<-merge(longdata1[,c("lon","lat","V1")],longdata2[,c("lon","lat","V2")],by=c("lon","lat"),all=T)
    temp_noV1<-temp[!is.na(temp$V1),]
    narows<-which(is.na(temp_noV1$V2))
    
    for (i in narows){
      temp_noV1[i,"V2"]<-interp_pt(temp_noV1[i,"lon"],temp_noV1[i,"lat"],temp,lon2,lat2)
    }
    longdata<-temp_noV1
  }else{
    longdata<-merge(longdata1[,c("lon","lat","V1")],longdata2[,c("lon","lat","V2")],by=c("lon","lat"))
  }

  #Create a cosine latitude column
  longdata$coslat<-cos(longdata$lat*pi/180)
  longdata$cosV1<-longdata$V1*longdata$coslat
  longdata$cosV2<-longdata$V2*longdata$coslat
  
  Vproduct<-longdata$cosV1*longdata$cosV2
  
  V1bar<-mean(longdata$cosV1)
  V2bar<-mean(longdata$cosV2)
  
  V1diff<-longdata$cosV1-V1bar
  V2diff<-longdata$cosV2-V2bar
  VdiffProduct<-V1diff*V2diff
  
  if (centered==TRUE){
    r<-sum(VdiffProduct)/sqrt(sum(V1diff*V1diff)*sum(V2diff*V2diff))    
  }else{
    r<-sum(Vproduct)/sqrt(sum(longdata$cosV1*longdata$cosV1)*sum(longdata$cosV2*longdata$cosV2))
  }

  return(r)
}