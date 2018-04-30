source("~/tempestextremes/test/sector_funcs.R")

args=commandArgs(trailingOnly=TRUE)
if (length(args)<3){
  stop("Argument 1: input directory Argument 2: output directory  Argument 3: dataset")
}
f_dir=args[1]
data=args[3]
sector=args[4]
season=args[5]
#for (data in c("ERA","MERRA")){
#f_dir=sprintf("~/block_r_data/new_data/%s_data/stats/",data)

#for (sector in c("NA","NC","NP","SA","SI","SP")){
#for (sector in c("NA")){
#  for (season in c("MAM","JJA","SON","DJF")){
    search_patt=sprintf("%s_*_%s_%s_*stats.txt",data,season,sector)
    files_masterlist=list.files(path=f_dir,pattern=glob2rx(search_patt))
    files_masterlist=paste(f_dir,files_masterlist,sep="/")
    nFiles=length(files_masterlist)
    
df_tot_stitch<-data.frame(var=character(),bnum=numeric(),tstep=numeric(),date=character(),hr=numeric(),
                          minlat=numeric(),maxlat=numeric(),
                          minlon=numeric(),maxlon=numeric(),
                          centlat=numeric(),centlon=numeric(),area=numeric(),
                          season=character(),sector=character())
bnum=0

    for (a in files_masterlist){
      #gather the info from filename
      print(sprintf("Opening file %s",a))
      dirsplit<-unlist(strsplit(a,split="/"))
      finfo<-unlist(strsplit(dirsplit[length(dirsplit)],split="_"))
      vartype=finfo[5]
      ynum<-finfo[2]
      mnum<-ifelse(season=="MAM","03",ifelse(
        season=="JJA","06",ifelse(
          season=="SON","09","12")
        )
      )
      orig<-sprintf("%s-%s-01",ynum,mnum)
      lines<-readLines(a)
      print(sprintf("There are %d lines in this file, type %s.",length(lines),vartype))
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
          t<-as.Date(nhours/24,origin=orig)
          #t<-as.Date(nhours/24,origin="1986-12-01")
          #print(t)
          hr<-nhours%%24
#          infoline[8]<-infoline[8]*6371000
          df_new<-data.frame(var=vartype,bnum=bnum,tstep=infoline[1],date=t,hr=hr,
                             minlat=infoline[2],maxlat=infoline[3],
                             minlon=infoline[4],maxlon=infoline[5],
                             centlat=infoline[6],centlon=infoline[7],area=infoline[8],
                             season=season,sector=sector)

          df_tot_stitch<-rbind(df_tot_stitch,df_new)
        }
      }
      print(sprintf("df is now %d long",nrow(df_tot_stitch)))
    }
    print(sprintf("There were %s blobs in %s %s",sector,season,bnum))

if (data=="ERA"){
  df_tot_stitch$minlon_c<-lon_convert(df_tot_stitch$minlon)
  df_tot_stitch$maxlon_c<-lon_convert(df_tot_stitch$maxlon)
  df_tot_stitch$centlon_c<-lon_convert(df_tot_stitch$centlon)
}
if (data=="MERRA"){
  df_tot_stitch$minlon_c<-df_tot_stitch$minlon
  df_tot_stitch$maxlon_c<-df_tot_stitch$maxlon
  df_tot_stitch$centlon_c<-df_tot_stitch$centlon
  df_tot_stitch$minlon<-lon_convert2(df_tot_stitch$minlon_c)
  df_tot_stitch$maxlon<-lon_convert2(df_tot_stitch$maxlon_c)
  df_tot_stitch$centlon<-lon_convert2(df_tot_stitch$centlon_c)
}
#Save the dataframe
df_tot_stitch<-df_tot_stitch[format(df_tot_stitch$date,"%m_%d")!="02_29",]
statsname<-sprintf("%s/%s_%s_%s_stats_stitch_table.RData",args[2],data,season,sector)
save(list=c("df_tot_stitch"),file=statsname)
#  }
#}
#}
