source("~/tempestextremes/test/sector_funcs.R")



f_dir="/Volumes/ExFAT_drive/ERA_files/ERA_blobs/"
#search_patt="ERA_1986_DJF_NA_*stats.txt"
search_patt="ERA_2003_JJA_NA_*stats.txt"
files_masterlist=list.files(path=f_dir,pattern=glob2rx(search_patt))
files_masterlist=paste(f_dir,files_masterlist,sep="")
nFiles=length(files_masterlist)

df_tot<-data.frame(var=character(),bnum=numeric(),tstep=numeric(),date=character(),hr=numeric(),
                   minlat=numeric(),maxlat=numeric(),
                   minlon=numeric(),maxlon=numeric(),
                   centlat=numeric(),centlon=numeric(),area=numeric())
bnum=0
for (a in files_masterlist){
  #gather the info from filename
  finfo<-unlist(strsplit(a,split="_"))
  vartype=ifelse(finfo[length(finfo)]=="stats.txt","PV",ifelse(
    finfo[length(finfo)]=="Zstats.txt","Z","GHG"
  ))
  lines<-readLines(a)
  for (l in lines){
    #print(l)
    if (!is.na(pmatch("Blob",l))){
      #Split blob line
      blobline=unlist(strsplit(l,split=" "))
      #bnum=as.numeric(blobline[2])
      bnum=bnum+1
    }else{
      #read in info
      #time step,minlat,maxlat,minlon,maxlon,centlat,centlon,area
      infoline<-as.numeric(unlist(strsplit(l,split="\t")))
      nhours<-infoline[1]*6
      #print(sprintf("bnum is %f, nhours is %f,mult by 6 is %f.",bnum,infoline[1],nhours))
      t<-as.Date(nhours/24,origin="2003-06-01")
      #print(t)
      hr<-nhours%%24
      infoline[8]<-infoline[8]*6371000
      df_new<-data.frame(var=vartype,bnum=bnum,tstep=infoline[1],date=t,hr=hr,
                         minlat=infoline[2],maxlat=infoline[3],
                         minlon=infoline[4],maxlon=infoline[5],
                         centlat=infoline[6],centlon=infoline[7],area=infoline[8])
      df_tot<-rbind(df_tot,df_new)
    }
  }
}
df_tot$minlon_c<-lon_convert(df_tot$minlon)
df_tot$maxlon_c<-lon_convert(df_tot$maxlon)
df_tot$centlon_c<-lon_convert(df_tot$centlon)


