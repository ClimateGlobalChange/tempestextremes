lon_convert<-function(lon){
  distFrom180=lon-180.
  return(ifelse(
    distFrom180<0,
    lon,
    -(180-distFrom180)
  ))
}


#-180 to 180 -> 
lon_convert2<-function(lon){
  return(ifelse(lon<0,360+lon,lon))
}


merge_dfs<-function(df_stitch,df_nostitch,rfn="",textfn="",csvfn="",df_merged_name="",
                    byvec=c("datehour","area","var")){
  df_merged_name<-ifelse(df_merged_name=="","df_merged",df_merged_name)
	df_names<-names(df_stitch)
	#This data frame only has common rows
	df_comm<-merge(df_stitch,df_nostitch,by=byvec)
	#This data frame has all rows, both merged and not
	df_tot<-merge(df_stitch,df_nostitch,by=byvec,all=T)
	#This data frame has all the rows of the merged blobs
	df_istot<-df_tot[is.na(df_tot$bnum.y),]
	#This data frame has all the rows of the multiple blobs that make up the merged blob
	df_isnot<-df_tot[is.na(df_tot$bnum.x),]
	
	df_comm$bnum2<-df_comm$bnum.x
	df_tot$bnum2<-df_tot$bnum.x
	df_istot$bnum2<-df_istot$bnum.x
	df_isnot$bnum2<-df_isnot$bnum.y
	
	#Now checking blobs where there is not a match:
	#Check whether or not blobs from DetectBlobs are occurring at the same time step 
	#within the extent of the original StitchBlobs output
	for (t in unique(df_istot$datehour)){
		df_check<-df_istot[df_istot$datehour==t,]
		df_otherblobs<-df_isnot[df_isnot$datehour==t,]
		if (nrow(df_otherblobs)>0){
			for (n in 1:nrow(df_check)){
				#Get the min/max lat and lon extent
				#Might need to deal with periodic boundary condition!
				clatmin<-df_check[n,"minlat.x"]
				clatmax<-df_check[n,"maxlat.x"]
				PER_BOUND<-FALSE
				clonmin<-df_check[n,"minlon.x"]
				clonmax<-df_check[n,"maxlon.x"]
			  axis180<-ifelse((clonmin<0 | clonmax<0),TRUE,FALSE)
				if (clonmin>clonmax){
					PER_BOUND<-TRUE
					clonmin<-ifelse(axis180==FALSE,lon_convert(clonmin),lon_convert2(clonmin))
					clonmax<-ifelse(axis180==FALSE,lon_convert(clonmax),lon_convert2(clonmax))
				}
				for (y in 1:nrow(df_otherblobs)){
					#Does it fall within the bounds of the big stitched blob?
					blatmin<-df_otherblobs[y,"minlat.y"]
					blatmax<-df_otherblobs[y,"maxlat.y"]
					blonmin<-df_otherblobs[y,"minlon.y"]
					blonmax<-df_otherblobs[y,"maxlon.y"]
					if (PER_BOUND==TRUE){
					  blonmin<-ifelse(axis180==FALSE,lon_convert(blonmin),lon_convert2(blonmin))
					  blonmax<-ifelse(axis180==FALSE,lon_convert(blonmax),lon_convert2(blonmax))
					}
					if (((blatmin>=clatmin)&(blatmax<=clatmax)&
						(blonmin>=clonmin)&(blonmax<=clonmax)&
						(df_check[n,"var"]==df_otherblobs[y,"var"]))){
							#Replace the unmerged blob number with the merged blob number
							df_otherblobs[y,"bnum.x"]<-df_check[n,"bnum.x"]
							df_comm<-rbind(df_comm,df_otherblobs[y,])				
					}
				}	
			}	
		}
	}
	#Clean up the output
	tcol<-grep("datehour",names(df_comm))
	vcol<-grep("var",names(df_comm))
	bcol<-grep("bnum.x",names(df_comm))
	b2col<-grep("bnum2",names(df_comm))
	acol<-grep("area",names(df_comm))
	akcol<-grep("area_km",names(df_comm))
	acol<-acol[!acol %in% akcol]
	#Get all of the columns with .y in the name
	ycol<-grep(".y",names(df_comm))
	bycol<-grep("bnum.y",names(df_comm))
	#Remove bnum.y from columns
	ycol<-ycol[!ycol %in% bycol]

	
	df_return<-df_comm[,c(tcol,ycol,bcol,b2col,acol,vcol)]
	colnames(df_return)<-gsub("\\.y","",names(df_return))
	colnames(df_return)<-gsub("\\.x","",names(df_return))
	df_final<-df_return[,c(df_names,"bnum2")]
	#Switch file and bnum2
	nlast<-length(names(df_final))
	df_final<-df_final[,c(1:(nlast-2),nlast,nlast-1)]
	df_final<-df_final[order(df_final$datehour),]
	if (rfn!=""){
	  assign(df_merged_name,df_final)
	  assign("df_name",df_merged_name)
	  save(list=c(df_merged_name,"df_name"),file=rfn) 
	}
	if (textfn!=""){
	  write.table(df_final,file=textfn,sep="\t",row.names=FALSE,quote=FALSE)
	}
	if (csvfn!=""){
	  write.csv(df_final,file=csvfn,row.names=FALSE,quote=FALSE)
	}
	#Return the merged data frame
	return(df_final)
}