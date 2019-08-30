combine_dfs<-function(flist,ftype,rfn="",textfn="",csvfn="",df_outname=""){
  df_data<-data.frame(NULL)
  df_outname<-ifelse(df_outname=="","df_data",df_outname)
  for (f in flist){
    if (ftype=="R"){
      load(f)
      df_tot<-get(df_name)
    }else if (ftype=="text"){
      df_tot<-read.table(f,header=TRUE,sep="\t")
    }else if (ftype=="CSV"){
      df_tot<-read.csv(f)
    }else{
      stop("Invalid file type.")
    }
    df_data<-rbind(df_data,df_tot)
  }
  if (rfn!=""){
    assign(df_outname,df_data)
    assign("df_name",df_outname)
    save(list=c(df_name,"df_name"),file=rfn) 
    print(sprintf("Wrote file %s",rfn))
  }
  if (textfn!=""){
    write.table(df_data,file=textfn,sep="\t",row.names=FALSE,quote=FALSE)
    print(sprintf("Wrote file %s",textfn))
  }
  if (csvfn!=""){
    write.csv(df_tot,file=csvfn,row.names=FALSE,quote=FALSE)
    print(sprintf("Wrote file %s",csvfn))
  }
  return(df_data)
}