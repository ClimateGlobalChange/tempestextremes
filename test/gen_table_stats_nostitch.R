source("~/tempestextremes/test/sector_funcs.R")

f_dir="/Volumes/ExFAT_drive/ERA_files/ERA_detect/"


for (sector in c("NA","NC","NP","SA","SI","SP")){
#for (sector in c("NA")){
  for (season in c("MAM","JJA","SON","DJF")){
    search_patt=sprintf("ERA_*_%s_%s_*stats_nostitch.txt",season,sector)
    files_masterlist=list.files(path=f_dir,pattern=glob2rx(search_patt))
    files_masterlist=paste(f_dir,files_masterlist,sep="")
    nFiles=length(files_masterlist)
    
    df_tot_nostitch<-data.frame(var=character(),bnum=numeric(),tstep=numeric(),date=character(),hr=numeric(),
                       minlat=numeric(),maxlat=numeric(),
                       minlon=numeric(),maxlon=numeric(),
                       centlat=numeric(),centlon=numeric(),area=numeric())
    bnum=0
    for (a in files_masterlist){
      #gather the info from filename
      print(sprintf("Opening file %s",a))
      finfo<-unlist(strsplit(a,split="_"))
      vartype=ifelse(finfo[length(finfo)-1]=="stats","PV",ifelse(
        finfo[length(finfo)-1]=="Zstats","Z","GHG"
      ))
      ynum<-finfo[length(finfo)-4]
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
          infoline[8]<-infoline[8]*6371000
          df_new<-data.frame(var=vartype,bnum=bnum,tstep=infoline[1],date=t,hr=hr,
                             minlat=infoline[2],maxlat=infoline[3],
                             minlon=infoline[4],maxlon=infoline[5],
                             centlat=infoline[6],centlon=infoline[7],area=infoline[8])
          df_tot_nostitch<-rbind(df_tot_nostitch,df_new)
        }
      }
      print(sprintf("df is now %d long",nrow(df_tot_nostitch)))
    }
    print(sprintf("There were %s blobs in %s %s",sector,season,bnum))
    df_tot_nostitch$minlon_c<-lon_convert(df_tot_nostitch$minlon)
    df_tot_nostitch$maxlon_c<-lon_convert(df_tot_nostitch$maxlon)
    df_tot_nostitch$centlon_c<-lon_convert(df_tot_nostitch$centlon)
    #Save the dataframe
    df_tot_nostitch<-df_tot_nostitch[format(df_tot_nostitch$date,"%m_%d")!="02_29",]
    statsname<-sprintf("~/block_r_data/stats_nostitch_%s_%s_table.RData",season,sector)
    save(list=c("df_tot_nostitch"),file=statsname)
  }
}

