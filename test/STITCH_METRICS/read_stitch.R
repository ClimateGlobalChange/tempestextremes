# This function reads the BlobStats files into a single data table
read_stats_to_table<-function(flist,nhrs,var="",
                              rfn="",textfn="",csvfn=""){
  nres<-nhrs/24
  fnum<-1
  df_tot<-data.frame(NULL)
  nline<-1
  for (f in flist){
    #Open each BlobStats file
    #The following is predicated upon Paul changing the BlobStats output!!
    #Might need to change the first few lines based on what the resulting changes are
    
    #First line: The start date at time 0 (in format YYYY-MM-DD-HH)
    fl<-readLines(f)
    startdate<-unlist(strsplit(fl[1],split="-"))
    orig<-substr(fl[1],1,10)
    #Currently assuming that nres is manually entered
    #Second line: Column names
    varnames<-unlist(strsplit(fl[2],split="\t"))
    tname<-varnames[1]
    
    fdat<-fl[3:length(fl)]
    for (l in fdat){
      if (!is.na(pmatch("Blob",l))){
        #Split the blobline
        blobline<-unlist(strsplit(l,split=" "))
        bnum<-blobline[2]
      }else{
        #Read in the info
        df_line<-data.frame(matrix(NA,nrow=1,ncol=length(varnames)))
        colnames(df_line)<-varnames
        infoline<-as.numeric(unlist(strsplit(l,split="\t")))
        #print(varnames)
        #print(length(varnames))
        #print(infoline)
        #print(length(infoline))
        #print(infoline)
        df_line[1,varnames]<-infoline
        #print(df_line)
        #print(df_line[1,tname])

        date<-as.Date(as.numeric(df_line[1,tname])*nres,origin=orig)
        hr<-(df_line[1,tname]*nhrs)%%24
        df_tot[nline,"datehour"]<-sprintf("%s %02d:00:00",as.character(date),hr)
        for (v in 2:length(varnames)){
          df_tot[nline,varnames[v]]<-as.numeric(df_line[1,v])
        }
        if (!is.null(df_line[1,"area"])){
          df_tot[nline,"area_km"]<-df_line[1,"area"]*4*pi*(6371^2)
        }
        if (var!=""){
          df_tot[nline,"var"]<-var
        }
        df_tot[nline,"bnum"]<-bnum
        df_tot[nline,"file"]<-f
        nline<-nline+1
      }
    }
    fnum<-fnum+1
  }

  if (rfn!=""){
    save(list=c("df_tot"),file=rfn) 
  }
  if (textfn!=""){
    write.table(df_tot,file=textfn,sep="\t",row.names=FALSE,quote=FALSE)
  }
  if (csvfn!=""){
    write.csv(df_tot,file=csvfn,row.names=FALSE,quote=FALSE)
  }
  return(df_tot)
}