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
#Interpolates point either along x, y, or bilinearly
#Note that it ALWAYS interpolates from 2's grid to 1's grid
# interp_pt<-function(x, y, df, lon2, lat2){
#   #Grab the 4 points and their values from the lon2 and lat2 axis
#   #Find the closest x value on the lon2 axis
#   #print(sprintf("x: %f y: %f",x,y))
#   xdist<-sqrt((x-lon2)^2)
#   xi<-which(xdist==min(xdist))
#   #print(sprintf("xi: %d",xi))
#   if (length(xi)>1){
#     #These will be the left and right points
#     xl<-round(lon2[xi[1]],3)
#     xr<-round(lon2[xi[2]],3)
#   }else{
#     #Is this point to the right or left of x?
#     #It's to the right
#     if (lon2[xi]>x){
#       if (xi==1){
#         xr<-round(x,3)
#         xl<-round(x,3)
#       }
#       else{
#         xr<-round(lon2[xi],3)
#         xl<-round(lon2[xi-1],3)
#       }
#       #it's to the left
#     }else if (lon2[xi]<x){
#       xl<-round(lon2[xi],3)
#       xr<-round(lon2[xi+1],3)
#     }else if(lon2[xi]==x){
#       xl<-round(x,3)
#       xr<-round(x,3)
#     } 
#   }
#   #print(sprintf("xl,xr: %f,%f",xl,xr))
#   #Now the y points
#   ydist<-sqrt((y-lat2)^2)
#   yi<-which(ydist==min(ydist))
#   #print(sprintf("yi: %d",yi))
#   if (length(yi)>1){
#     yb<-round(lat2[yi[1]],3)
#     yt<-round(lat2[yi[2]],3)
#   }else{
#     #It's above the y point
#     if (lat2[yi]>y){
#       yt<-round(lat2[yi],3)
#       yb<-round(lat2[yi-1],3)
#     }else if (lat2[yi]<y){
#       yb<-round(lat2[yi],3)
#       yt<-round(lat2[yi+1],3)
#     }else if (lat2[yi]==y){
#       yb<-round(y,3)
#       yt<-round(y,3)
#     }
#   }
#   #print(sprintf("yt,yb: %f,%f",yt,yb))
#   #Get the distance from the interpolated point and average
#   ipt<-0
#   if (xl==xr && yt==yb){
#     #print("same point")
#     ipt<-df[df$lon==xl & df$lat==yt,"V1"]
#   }else if (xl==xr){
#     #Only interpolate in the y direction
#     #print("Interpolating y only")
#     ipt<-((yt-y)/(yt-yb))*df[df$lon==x & df$lat==yb,"V2"] + ((y-yb)/(yt-yb))*df[df$lon==x & df$lat==yt,"V2"]
#   }else if (yt==yb){
#     #Only interpolate in the x direction
#     #print("Interpolating x only")
#     ipt<-((xr-x)/(xr-xl))*df[df$lon==xl & df$lat==y,"V2"] + ((x-xl)/(xr-xl))*df[df$lon==xr & df$lat==y,"V2"]
#   }else{
#     #print("Interpolating both x and y")
#     p11<-df[df$lon==xl & df$lat==yb,"V2"]
#     p12<-df[df$lon==xl & df$lat==yt,"V2"]
#     p21<-df[df$lon==xr & df$lat==yb,"V2"]
#     p22<-df[df$lon==xr & df$lat==yt,"V2"]
#     ipt<-(1/((xr-xl)*(yt-yb)))*(p11*(xr-x)*(yt-y) + p21*(x-xl)*(yt-y) + p12*(xr-x)*(y-yt) + p22*(x-xl)*(y-yb))
#   }
#   return(ipt)
# }

linint<-function(from_arr,from_lat,from_lon,to_lat,to_lon){
  to_lon<-round(to_lon,3)
  to_lat<-round(to_lat,3)
  from_lon<-round(from_lon,3)
  from_lat<-round(from_lat,3)
  dx<-to_lon[2]-to_lon[1]
  dy<-to_lat[2]-to_lat[1]
  print(sprintf("%f %f",dx,dy))
  new_grid<-bilinear.grid(x=from_lon,y=from_lat,z=from_arr,
                          xlim=range(to_lon),ylim=range(to_lat),
                          dx=dx,dy=dy)
  #print(dim(new_grid$z))
  #print(dim(from_arr))
  return(new_grid$z)
}

#Finds the spatial similarity between two blobs
similarity_weighted_old<-function(arr1,arr2,lats,lons,lats2=NULL,lons2=NULL,regrid=FALSE){
  longdata1<-melt(arr1,value.name = "V1")
  longdata1$V1[longdata1$V1>0]<-1
	longdata1$lon<-round(lons[longdata1$Var1],3)
  longdata1$lat<-round(lats[longdata1$Var2],3)
 # print("Melting arr2")
  longdata2<-melt(arr2,value.name="V2")
  longdata2$V2[longdata2$V2>0]<-1
	if (is.null(lons2)){
		longdata2$lon<-round(lons[longdata2$Var1],3)
	}else{
		longdata2$lon<-round(lons2[longdata2$Var1],3)
	}
	if (is.null(lats2)){
		longdata2$lat<-round(lats[longdata2$Var2],3)
	}else{
		longdata2$lat<-round(lats2[longdata2$Var2],3)
	}
  #NOTE THAT REGRID WILL ALWAYS BE FROM 2 TO 1
  if (regrid==TRUE){
    #Find lat/lon missing pairs
    #Basically wherever V2 is NA
    temp<-merge(longdata1[,c("lon","lat","V1")],longdata2[,c("lon","lat","V2")],by=c("lon","lat"),all=T)
    #Omit all instances where the coarser grid is NA
    temp_noV1<-temp[!is.na(temp$V1),]
    #Get all the row numbers where there's an NA
    narows<-which(is.na(temp_noV1$V2))
    # for (i in narows){
    #   temp_noV1[i,"V2"]<-interp_pt(temp_noV1[i,"lon"],temp_noV1[i,"lat"],temp,lons2,lats2)
    # }
    temp_noV1[narows,"V2"]<-sapply(narows,function(i){
      interp_pt(temp_noV1[i,"lon"],temp_noV1[i,"lat"],temp,lons2,lats2)
    })
    temp_noV1$V2<-ifelse(temp_noV1$V2>0.3,1,0)
    longdata<-temp_noV1
  }else{
    longdata<-merge(longdata1[,c("lon","lat","V1")],longdata2[,c("lon","lat","V2")],by=c("lon","lat"))
  }
  #print(longdata)

  longdata$V1_w<-longdata$V1*cos(longdata$lat*pi/180)
  longdata$V2_w<-longdata$V2*cos(longdata$lat*pi/180)
  intersect_12<-longdata[longdata$V1>0 & longdata$V2>0,]
  ncommon<-nrow(intersect_12)
  
  union_12<-longdata[longdata$V1>0 | longdata$V2>0,]
  union_12$VU_w<-cos(union_12$lat*pi/180)
  nunion<-nrow(union_12)
  #print(sprintf("There are %d common points for union (%d points)",ncommon,nunion))
  
  sum_iw<-sum(intersect_12$V1_w)
  sum_uw<-sum(union_12$VU_w)
  print(sprintf("Weighted sums are %f and %f",sum_iw,sum_uw))
  if (abs(sum_uw)<0.00000000001){
    ratio<-0
  }else{
    ratio<-sum_iw/sum_uw
  }
  print(sprintf("ratio is %f",ratio))
  return(ratio)
}

similarity_weighted<-function(arr1,arr2,lats,lons,lats2=NULL,lons2=NULL,regrid=FALSE){
  longdata1<-melt(arr1,value.name = "V1")
  longdata1$V1<-ifelse(longdata1$V1>0,1,0)
  longdata1$lon<-round(lons[longdata1$Var1],3)
  longdata1$lat<-round(lats[longdata1$Var2],3)
  
  if (is.null(lats2)){
    lats2<-lats
  }
  if (is.null(lons2)){
    lons2<-lons
  }
 
  if (regrid==TRUE){
    arr2_analyze<-linint(arr2,lats2,lons2,lats,lons)
  }else{
    arr2_analyze<-arr2
  }
  
  #print(dim(arr1))
  #print(dim(arr2_analyze))
  longdata2<-melt(arr2_analyze,value.name="V2")
  longdata2$V2<-ifelse(longdata2$V2>0.3,1,0)
  longdata2$lon<-round(lons[longdata2$Var1],3)
  longdata2$lat<-round(lats[longdata2$Var2],3)
  longdata<-merge(longdata1,longdata2,by=c("lon","lat"))
  
  
  longdata$V1_w<-longdata$V1*cos(longdata$lat*pi/180)
  longdata$V2_w<-longdata$V2*cos(longdata$lat*pi/180)
  intersect_12<-longdata[longdata$V1>0 & longdata$V2>0,]
  ncommon<-nrow(intersect_12)
  
  union_12<-longdata[longdata$V1>0 | longdata$V2>0,]
  union_12$VU_w<-cos(union_12$lat*pi/180)
  nunion<-nrow(union_12)
  #print(sprintf("There are %d common points for union (%d points)",ncommon,nunion))
  
  sum_iw<-sum(intersect_12$V1_w)
  sum_uw<-sum(union_12$VU_w)
  print(sprintf("Weighted sums are %f and %f",sum_iw,sum_uw))
  if (abs(sum_uw)<0.00000000001){
    ratio<-0
  }else{
    ratio<-sum_iw/sum_uw
  }
  print(sprintf("ratio is %f",ratio))
  return(ratio)
}

#Finds the probability of overlap using the original data frame and the data frame of overlaps
prob_calc<-function(df,dfo,V1,V2){
  V1count<-nrow(df[df$var==V1,])
  V2count<-nrow(df[df$var==V2,])
  #print(sprintf("V1: %d, V2:%d",V1count,V2count))
  V12count<-0
  for (d in sort(unique(dfo$datehour))){
    dsub<-dfo[dfo$datehour==d,]
    #print(dsub)
    nV1<-length(unique(dsub$V1bnum2))
    nV2<-length(unique(dsub$V2bnum2))
    nV12<-nrow(dsub)
    #print(sprintf("V1: %d V2: %d, V12: %d",nV1,nV2,nV12))
    V12count<-V12count+min(nV1,nV2,nV12)
  }
  p1given2<-V12count/V2count
  #print(p1given2)
  p2given1<-V12count/V1count
  #print(p2given1)
  return(c(p1given2,p2given1))
}
#Finds the Pearson pattern correlation between the averaged climatologies of the two datasets
pearson_arr<-function(arr1,arr2,lat1,lat2,lon1,lon2,interp=FALSE,centered=FALSE){
  
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

#Calculates the rmse of 2 wrt 1
rmse_calc<-function(arr1,arr2,lat1,lat2,lon1,lon2,interp=FALSE){
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
  #print(length(diff_data2))
  #rint(nrow(longdata))
  #print(range(diff_data2))
  ame<-sum(abs(longdata$cosV1-longdata$cosV2))/nrow(longdata)
  rmse<-sqrt(sum(diff_data2)/length(diff_data2))
  print(sprintf("AME is %f and RMSE is %f",ame,rmse))
  return(rmse)
}

#Finds the index of closest value in a vector
nearest_ind<-function(val,vec){
  dist_vec<-abs(val-vec)
  i<-which(dist_vec==min(dist_vec))
  return(i[1])
}

#This function looks at overlaps between blobs that came from different StitchBlobs outputs
#at the same time step. Requires 2 distinct variable names in the var column for the df
# Will probably require RData file for the similarity calculations, too much of a pain
# in the behind otherwise
overlaps_calc<-function(df1,blobs1,time_format1,lat1,lon1,blobs2,
                        df2=NULL,time_format2=NULL,lat2=NULL,lon2=NULL,
                        rfn_ps="",txt_overlaps="",txt_ps="",regrid=FALSE){

	df_overlaps<-data.frame(NULL)
	#Default assumption that all of the data is contained within a single data frame
	df_analyze<-df1
	#If not, combine data frames
	if (!is.null(df2)){
		df_analyze<-rbind(df_analyze,df2)
	}
	
	#If not adding a second set of values for time and spatial resolution,
	# make them non-null
	if (is.null(time_format2)){
		time_format2<-time_format1
	}
	if (is.null(lat2)){
		lat2<-lat1
	}
	if (is.null(lon2)){
		lon2<-lon1
	}
	
	#Does the longitude axis have the same range of values?
	#Create longitude axes that have both 0->360 and -180->180 axes
	lon1_180<-lon_convert(lon1)
	lon2_180<-lon_convert(lon2)
	lon1_360<-lon_convert2(lon1)
	lon2_360<-lon_convert2(lon2)
	
	df_analyze$minlon_180<-lon_convert(df_analyze$minlon)
	df_analyze$minlon_360<-lon_convert2(df_analyze$minlon)
	df_analyze$maxlon_180<-lon_convert(df_analyze$maxlon)
	df_analyze$maxlon_360<-lon_convert2(df_analyze$maxlon)
	df_analyze$centlon_180<-lon_convert(df_analyze$centlon)
	df_analyze$centlon_360<-lon_convert2(df_analyze$centlon)

	#Use 180 or 360?
	if (lon1_360[1]>lon1_360[length(lon1_360)] | lon2_360[1]>lon2_360[length(lon2_360)]){
		lon1_analyze<-lon1_180
		lon2_analyze<-lon2_180
		minlon_analyze<-"minlon_180"
		maxlon_analyze<-"maxlon_180"
		centlon_analyze<-"centlon_180"
	}else{
		lon1_analyze<-round(lon1_360,3)
		lon2_analyze<-round(lon2_360,3)
		minlon_analyze<-"minlon_360"
		maxlon_analyze<-"maxlon_360"
		centlon_analyze<-"centlon_360"
	}
	
	#Restrict df values to the range of the longitude axis
	lon_min<-min(lon1_analyze,lon2_analyze)
	lon_max<-max(lon1_analyze,lon2_analyze)
	
	df_analyze<-df_analyze[(df_analyze[,minlon_analyze]>=lon_min &
		 df_analyze[,maxlon_analyze]<=lon_max),]
	
	#If the bnum2 column doesn't exist, create a column with the bnum2 name
	if (! "bnum2" %in% names(df_analyze)){
		df_analyze$bnum2<-df_analyze$bnum
	}

	#Get the variable names of the two algorithms being compared
	varNames<-unique(as.character(df_analyze$var))
	if (length(varNames)<2){
		stop("There is only one value in var column. Check input data frame.")
	}
	V1<-varNames[1]
	V2<-varNames[2]
	
	print(sprintf("V1 is %s and V2 is %s",V1,V2))
	sim_time<-data.frame(NULL)
	#Take each unique time step in the data frame 
	#Find the common time steps between the df and the NetCDFs
	timeNetCDF<-intersect(time_format1,time_format2)
	time_intersect<-intersect(df_analyze$datehour,timeNetCDF)
	df_analyze<-df_analyze[df_analyze$datehour %in% time_intersect,]
	
	nr<-1
	
  #Get the time indices that are in the NetCDF intersect for V1
	time_sub1<-which(time_format1 %in% timeNetCDF)
	#Get the time indices that are in the NetCDF intersect for V2
	time_sub2<-which(time_format2 %in% timeNetCDF)
	
	#First, find the pearson correlation between the averaged patterns
	b1avg<-apply(blobs1[,,time_sub1],c(1,2),mean)
	b2avg<-apply(blobs2[,,time_sub2],c(1,2),mean)
	
	pearson_num<-pearson_arr(b1avg,b2avg,round(lat1,3),round(lat2,3),round(lon1_analyze,3),round(lon2_analyze,3),interp=regrid)
	rmse_num<-rmse_calc(b1avg,b2avg,round(lat1,3),round(lat2,3),round(lon1_analyze,3),round(lon2_analyze,3),interp=regrid)
	
	for (d in sort(unique(time_intersect))){
	  print(d)
		dsub<-df_analyze[df_analyze$datehour==d,]
		dsub<-dsub[!duplicated(dsub),]
    #print(dsub)
		#How many unique types of blob variables are there at the time?
		varname<-as.character(sort(unique(as.character(dsub$var))))
		nvar<-length(varname)
		if (nvar>1){
			df1<-dsub[dsub$var==V1,]
			df2<-dsub[dsub$var==V2,]
			d1i<-which(time_format1==d)
			#time1_inds<-c(time1_inds,d1i)
			d2i<-which(time_format2==d)
			#time2_inds<-c(time2_inds,d2i)
			#Loop through the variables
			for (b1 in 1:nrow(df1)){
				#Find the blob boundaries
			  toplat1<-df1[b1,"maxlat"]
			  botlat1<-df1[b1,"minlat"]
			  rlon1<-df1[b1,maxlon_analyze]
			  llon1<-df1[b1,minlon_analyze]
				it1<-nearest_ind(toplat1,lat1)
				ib1<-nearest_ind(botlat1,lat1)
				il1<-nearest_ind(rlon1,lon1_analyze)
				ir1<-nearest_ind(llon1,lon1_analyze)
				#print(sprintf("Bounds (index) for V1: %d %d %d %d",it1,ib1,il1,ir1))
				#Grab the time slice
				v1slice<-blobs1[,,d1i]
				#Zero out everything but the blob in the longitude direction
				v1slice[c(1:min(il1,ir1),max(il1,ir1):dim(blobs1)[1]),]<-0
				#Zero out everything but the blob in the latitude direction
				v1slice[,c(1:min(it1,ib1),max(it1,ib1):dim(blobs1)[2])]<-0
				for (b2 in 1:nrow(df2)){
				  #Do the boundaries even overlap?
				  toplat2<-df2[b2,"maxlat"]
				  botlat2<-df2[b2,"minlat"]
				  rlon2<-df2[b2,maxlon_analyze]
				  llon2<-df2[b2,minlon_analyze]
				  print(sprintf("bounds for 1: %f,%f,%f,%f",toplat1,botlat1,rlon1,llon1))
				  print(sprintf("bounds for 2: %f,%f,%f,%f",toplat2,botlat2,rlon2,llon2))
				  has_overlap<-TRUE
				  if (rlon1<llon2){
				    print("1 is left of 2")
				    has_overlap<-FALSE
				  }
				  if (rlon2<llon1){
				    print("1 is right of 2")
				    has_overlap<-FALSE
				  }
				  if (botlat1>toplat2){
				    print("1 is above 2")
				    has_overlap<-FALSE
				  }
				  if (botlat2>toplat1){
				    print("1 is below 2")
				    has_overlap<-FALSE
				  }
          if (has_overlap==TRUE){
            it2<-nearest_ind(toplat2,lat2)
            ib2<-nearest_ind(botlat2,lat2)
            il2<-nearest_ind(rlon2,lon2_analyze)
            ir2<-nearest_ind(llon2,lon2_analyze)
            #print(sprintf("Bounds (index) for V2: %d %d %d %d",it2,ib2,il2,ir2))
            v2slice<-blobs2[,,d2i]
            v2slice[c(1:min(il2,ir2),max(il2,ir2):dim(blobs2)[1]),]<-0
            v2slice[,c(1:min(it2,ib2),max(it2,ib2):dim(blobs2)[2])]<-0
            s12<-similarity_weighted(v1slice,v2slice,lat1,lon1_analyze,lat2,lon2_analyze,regrid)
            if (s12>0){
              df_overlaps[nr,"datehour"]<-d
              df_overlaps[nr,"similarity"]<-s12
              col_copy<-c("bnum","bnum2","minlat","maxlat",minlon_analyze,
                          maxlon_analyze,"centlat",centlon_analyze)
              col_names<-c("bnum","bnum2","minlat","maxlat","minlon",
                           "maxlon","centlat","centlon")
              for (v in 1:length(col_copy)){
                df_overlaps[nr,sprintf("V1%s",col_names[v])]<-df1[b1,col_copy[v]]
              }
              for (v in 1:length(col_copy)){
                df_overlaps[nr,sprintf("V2%s",col_names[v])]<-df2[b2,col_copy[v]]
              }
              nr<-nr+1
            }
          }

				}
			}
		}
	}
  probs<-prob_calc(df_analyze,df_overlaps,V1,V2)
  p1given2<-probs[1]
  p2given1<-probs[2]
  sim_25<-quantile(df_overlaps$similarity,0.25)
  sim_50<-quantile(df_overlaps$similarity,0.5)
  sim_75<-quantile(df_overlaps$similarity,0.75)
  
  #Find the averaged blocking climatologies of the two datasets
  
  
  saved_names<-c("V1","V2","p1given2","p2given1","sim_25","sim_50","sim_75",
                 "df_overlaps","pearson_num","rmse_num")
  return_list<-list()
  for (v in saved_names){
    return_list[[v]]<-get(v)
  }
  if (rfn_ps!=""){
    save(list=saved_names,file=rfn_ps)
  }
  if (txt_overlaps!=""){
    write.table(df_overlaps,file=txt_overlaps,sep="\t",row.names=FALSE,quote=FALSE)
  }
  if (txt_ps!=""){
    sink(txt_ps)
    cat(sprintf("Variable 1: %s \n",V1))
    cat(sprintf("Variable 2: %s \n",V2))
    cat(sprintf("Pearson correlation between average of 1 and 2: %f \n",pearson_num))
    cat(sprintf("RMSE between average of 1 and 2: %f \n",rmse_num))
    cat(sprintf("Probability of %s given %s: %f \n",V1,V2,p1given2))
    cat(sprintf("Probability of %s given %s: %f \n",V2,V1,p2given1))
    cat(sprintf("25th percentile similarity value: %f \n",sim_25))
    cat(sprintf("50th percentile similarity value: %f \n",sim_50))
    cat(sprintf("75th percentile similarity value: %f \n",sim_75))
    sink()
  }
  return(return_list)
}