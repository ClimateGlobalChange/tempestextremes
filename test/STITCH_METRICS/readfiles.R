# This function reads the BlobStats files into a single data table
read_stats_to_table<-function(flist,var="",
                              rfn="",textfn="",csvfn="",df_name=""){

  df_tot<-data.frame(NULL)
  df_name<-ifelse(df_name=="","df_tot",df_name)
  nline<-1
  for (f in flist){
    print(sprintf("Reading in %s",f))
    #Open each BlobStats file
    fl<-readLines(f)
    if (length(fl)>1){
      #First line: Column names
      varnames<-unlist(strsplit(fl[1],split=","))
      if (length(varnames)<2){
        stop(sprintf("Check that file %s has the correct headers",f))
      }
      tname<-varnames[1]
      varnames<-varnames[2:length(varnames)]
      fdat<-fl[2:length(fl)]
      for (l in fdat){
        if (!is.na(pmatch("Blob",l))){
          #Split the blobline
          blobline<-unlist(strsplit(l,split="\\s+"))
          bnum<-blobline[2]
        }else{
          #Read in the info
          df_line<-data.frame(matrix(NA,nrow=1,ncol=length(varnames)))
          colnames(df_line)<-varnames
          #Line contains the time, then all the other data
          infoline<-unlist(strsplit(l,split="\\s+"))
          df_line[1,]<-as.numeric(infoline[2:length(infoline)])
          #Format the date string
          date_vec<-unlist(strsplit(infoline[1],split="-"))
          #Year, Month, Day, Seconds
          d<-paste(date_vec[1:3],collapse="-")
          hr<-as.integer(as.numeric(date_vec[4])/(60*60))
          df_tot[nline,"datehour"]<-sprintf("%s %02d:00:00",d,hr)
          
          for (v in 1:length(varnames)){
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
    }else{
      print("No blobs found")
    }
  }

  if (rfn!=""){
    assign(df_name,df_tot)
    
    save(list=c(df_name,"df_name"),file=rfn) 
    print(sprintf("Wrote file %s",rfn))
  }
  if (textfn!=""){
    write.table(df_tot,file=textfn,sep="\t",row.names=FALSE,quote=FALSE)
    print(sprintf("Wrote file %s",textfn))
  }
  if (csvfn!=""){
    write.csv(df_tot,file=csvfn,row.names=FALSE,quote=FALSE)
    print(sprintf("Wrote file %s",csvfn))
  }
  return(df_tot)
}