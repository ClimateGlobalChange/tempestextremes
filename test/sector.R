# This file generates the dataframe that contains all of the BlobStats data
# It looks for the wildcard pattern "statlist" (list of blobstat files for
# each dataset) and for each file on that list, reads the stat data into 
# the dataframe
########################
# Dataframe columns:  #
#######################
# Dataset: name of the dataset 
# nYears: number of years per dataset (necessary for normalizing histograms)
# Season: DJF, MAM, JJA, SON
# year: calendar year (or model run year)
# Blob: number of the blob within the file
# time: number of time steps since beginning of file
# minlat/maxlat: minimum/maximum latitudinal extent of blob
# minlon/maxlon: minimum/maximum longitudinal extent of blob
# centlat/centlon: location of blob centers
# area: current size of blob
# Hemi: hemisphere (NH/SH)
# Sector: Pacific, Atlantic, Indian, etc
# minlon_c/maxlon_c/centlon_c: longitude recalculated to 180W->180E

source("~/tempestextremes/test/sector_funcs.R")
#setwd("~/data_netcdf")
setwd("~/new_stats/")

latlon_dist<-function(lat1,lon1,lat2,lon2){
  dlat=(lat2-lat1)*pi/180
  
  #check periodicity
  if (abs(lon2-lon1>180)){
    lon1=lon_convert(lon1)
    lon2=lon_convert(lon2)
  }
  dlon=(lon2-lon1)*pi/180
  rlat1=lat1*pi/180
  rlat2=lat2*pi/180
  a=(sin(dlat/2))^2 + cos(rlat1)*cos(rlat2)*(sin(dlon/2))^2
  c=2*atan2(sqrt(a),sqrt(1-a))
  d=6371*c
  return(d)
}

#Data we have:
## Original variables
## Blobs
## Density of blobs
## Stats: max/min/center lat/lon, area
newblob_num=1
df_tot = data.frame(Dataset=character(),nYears=numeric(),Season=character(),year=numeric(),
                    Blob=numeric(),time=numeric(),minlat=numeric(),
                    maxlat=numeric(),minlon=numeric(),maxlon=numeric(),
                    centlat=numeric(),centlon=numeric(),
                    area=numeric(),Hemi=character(),Sector=character(),
                    minlon_c=numeric(),maxlon_c=numeric(),centlon_c=numeric(),
                    dist=numeric(),split=character(),newnum=numeric())

df_each=data.frame(Dataset=character(),Season=character(),year=numeric(),
                   Blob=numeric(),Hemi=character(),Sector_tot=character(),
                   centlat_start=numeric(),centlon_start=numeric(),
                   centlat_end=numeric(),centlon_end=numeric(),
                   centlat_span=numeric(),centlon_span=numeric(),duration=numeric(),
                   min_area=numeric(),max_area=numeric(),mean_area=numeric(),dist=numeric(),
                   sector_start=character(),sector_end=character(),split=character(),newnum=numeric(),Sector=character())

files_masterlist=list.files(pattern="*statlist")
print("Beginning to read in data files.")
for (a in files_masterlist){
  listfiles=readLines(a)
  nFiles=length(listfiles)
  
  #Find out total number of lines, + header lines
  tot_flines=0
  tot_hlines=0
  for (x in 1:nFiles){
    f=listfiles[x]
    #print(f)
    lines=readLines(f)
    nlines=length(lines)
    
    tot_flines = tot_flines + nlines
    #Grab all instances of header line
    BlobCounts=0
    for (i in 1:nlines){
      if (!is.na(pmatch("Blob",lines[i]))){
        BlobCounts = BlobCounts + 1
      }
    }
    tot_hlines = tot_hlines + BlobCounts
  }
  nrows=tot_flines-tot_hlines
  
  #Initialize data frame
  col = c('Dataset','nYears','Season','year','Blob','time','minlat',
          'maxlat','minlon','maxlon','centlat','centlon','area',
          'Hemi','Sector',"minlon_c","maxlon_c","centlon_c","dist","newnum")
  ncols=length(col)
  df=data.frame(matrix(ncol=ncols, nrow=nrows))
  colnames(df)<-col
  
  start_df_index=1
  end_df_index=0
  
  #Initialize summary dataframe
  col_summ=c('Dataset','Season','year','Blob','Hemi','Sector_tot',
             'centlat_start','centlon_start',
             'centlat_end','centlon_end',
             'centlat_span','centlon_span','duration',
             'min_area','max_area','mean_area',"dist",
             "sector_start", "sector_end","split","newnum","Sector")
  ncols_summ=length(col_summ)
  df_summ=data.frame(matrix(ncol=ncols_summ,nrow=1))
  colnames(df_summ)<-col_summ

  #Read in data!
  for (x in 1:nFiles){
    #Read file
    f=listfiles[x]
    print(f)
    #Get number of lines
    lines=readLines(f)
    nlines=length(lines)
    #Find all instances of header lines
    BlobIndices=c()
    for (i in 1:nlines){
      if (!is.na(pmatch("Blob",lines[i]))){
        BlobIndices=append(BlobIndices,i)
      }
    }
    
    #Grab year and season from file name
    fstring=unlist(strsplit(f,split="/"))
    ystring=unlist(strsplit(fstring[length(fstring)],split="_"))
    if(ystring[1]=="ERA"){
      dataname=ystring[1]
      year=as.numeric(ystring[2])
      season=ystring[3]
      sector_string=ystring[4]
      tstep=6
    }else{
      season=ystring[1]
      year=as.numeric(ystring[2])
      sector_string=ystring[3]
      tstep=3
      if (length(ystring)>5){
        datapaste=paste(ystring[length(ystring)-1],ystring[length(ystring)],sep="_")
        datastring=unlist(strsplit(datapaste,split="[.]"))
      }else{
        datastring=unlist(strsplit(ystring[length(ystring)],split="[.]"))
      }
      dataname=datastring[1]
    }
    
    #split strings and convert to data frame
    for (n in 1:length(BlobIndices)){
      if (n<(length(BlobIndices))){
        i1=BlobIndices[n]
        i2=BlobIndices[n+1]-1
      }else{
        i1=BlobIndices[n]
        i2=length(lines)
      }
      sub=lines[i1:i2]
      #number of lines (minus header)
      lsub=length(sub)-1
      #Number of the blob
      title=unlist(strsplit(sub[1],split=" "))
      nBlob=as.numeric(title[2])
      
      #Vector with Blob number, year number, season
      dvec=rep(dataname,lsub)
      
      bvec=rep(nBlob,lsub)
      yvec=rep(year,lsub)
      svec=rep(season,lsub)
      nYearVec=rep(nFiles,lsub)
      end_df_index=end_df_index + lsub
  
      #Data
      dat=do.call(rbind,strsplit(sub[2:length(sub)],split="\t"))
      class(dat)<-"numeric"
      #final_dat=cbind(dvec,nYearVec,svec,yvec,bvec,dat)

      #Add variables to summary dataframe
      df_summ[1,c('Blob','Dataset','Season','year')]=c(nBlob,dataname,season,year)
      
      df[start_df_index:end_df_index,c("Blob","year","nYears")] = sapply(cbind(bvec,yvec,nYearVec),as.numeric)
      df[start_df_index:end_df_index,c("Dataset","Season")] = cbind(dvec,svec)
      df[start_df_index:end_df_index,6:13]=dat

      df[start_df_index:end_df_index,"Hemi"]=ifelse(df[start_df_index:end_df_index,"centlat"]>0,"NH","SH")
      df[start_df_index:end_df_index,"Sector"]=sector_df(df[start_df_index:end_df_index,])
      #converted longitude degrees (from 0->360 to 180W->180E)
      
      df[start_df_index:end_df_index,"minlon_c"]=lon_convert(df[start_df_index:end_df_index,"minlon"])
      df[start_df_index:end_df_index,"maxlon_c"]=lon_convert(df[start_df_index:end_df_index,"maxlon"])
      df[start_df_index:end_df_index,"centlon_c"]=lon_convert(df[start_df_index:end_df_index,"centlon"])
      
      #Calculate distance traveled from timestep to timestep
      df[start_df_index,"dist"]=0
      df[(start_df_index+1):end_df_index,"dist"]=sapply((start_df_index+1):end_df_index, 
         function(x) latlon_dist(df[x-1,"centlat"],df[x-1,"centlon"],df[x,"centlat"],df[x,"centlon"]))
      subset=df[start_df_index:end_df_index,]
      df[start_df_index:end_df_index,"split"]=rep(ifelse(any(subset$dist>2000),"Y","N"),lsub)
      df[start_df_index:end_df_index,"newnum"]=rep(newblob_num,lsub)
      #print("finished filling in df")
      #Fill in rest of summary set

      df_summ[1,"Hemi"]=subset[1,"Hemi"]
      sub_sectors=unique(subset$Sector)
      if (length(sub_sectors)==1){
        df_summ[1,"Sector_tot"]=sub_sectors[1]
      }else{
        df_summ[1,"Sector_tot"]=paste(sub_sectors,collapse=" ")
      }
      df_summ[1,"centlat_start"]=subset$centlat[1]
      df_summ[1,"centlat_end"]=subset$centlat[length(subset$centlat)]
      df_summ[1,"centlon_start"]=subset$centlon[1]
      df_summ[1,"centlon_end"]=subset$centlon[length(subset$centlon)]
      
      df_summ[1,"sector_start"]=subset$Sector[1]
      df_summ[1,"sector_end"]=subset$Sector[length(subset$Sector)]
      
      df_summ[1,"min_area"]=min(subset$area)*6371*6371
      df_summ[1,"max_area"]=max(subset$area)*6371*6371
      df_summ[1,"mean_area"]=mean(subset$area)*6371*6371
      
      df_summ[1,"centlat_span"]=max(subset$centlat)-min(subset$centlat)
      #Deal with periodicity
      #Check if the Atlantic basin is within the sectors!

      cent_diff=max(subset$centlon)-min(subset$centlon)
      centc_diff=max(subset$centlon_c)-min(subset$centlon_c)      
      df_summ[1,"centlon_span"]=min(cent_diff,centc_diff)

      t=subset$time
      df_summ[1,"duration"]=tstep*(t[length(t)]-t[1])
      df_summ[1,"dist"]=sum(subset$dist)
      df_summ[1,"split"]=ifelse(any(subset$dist>2000),"Y","N")
      df_summ[1,"newnum"]=newblob_num
      df_summ[1,"Sector"]=sector_string
      df_each<-rbind(df_each,df_summ)
      
      start_df_index = end_df_index + 1
      newblob_num = newblob_num+1
    }
  }
  #df[,4:13]=lapply(df[4:13],as.numeric)
  #df$Hemi=ifelse(df$centlat>0,"NH","SH")
  #df$Sector=sector_df(df)
  

  #converted longitude degrees (from 0->360 to 180W->180E)
  #df$minlon_c=lon_convert(df$minlon)
  #df$maxlon_c=lon_convert(df$maxlon)
  #df$centlon_c=lon_convert(df$centlon)
  print(paste("finished with file",f))
  df_tot<-rbind(df_tot,df)
}

df_tot[,2]=lapply(df_tot[2],as.numeric)
df_tot$area=df_tot$area*6371*6371
df_each$days=df_each$duration/24.

colvec=c("blue","red","green","purple","pink")

df_each$col_data=ifelse(df_each$Dataset=="ERA","blueviolet",
                   ifelse(df_each$Dataset=="climo","blue",
                          ifelse(df_each$Dataset=="2xCO2","green",
                                 ifelse(df_each$Dataset=="SSTplus2","aquamarine2","cyan4"))))
df_each$short_name=ifelse(df_each$Dataset=="ERA","E",
                   ifelse(df_each$Dataset=="climo","CL",
                          ifelse(df_each$Dataset=="2xCO2","2C",
                                 ifelse(df_each$Dataset=="SSTplus2","S2","SC"))))


#test_vals$col_data=ifelse(test_vals$class=="ERA","orange",
#                        ifelse(test_vals$class=="climo","blue",
#                               ifelse(test_vals$class=="2xCO2","green",
#                                      ifelse(test_vals$class=="SSTplus2","pink","cyan4"))))

print("Finished generating data frame.")