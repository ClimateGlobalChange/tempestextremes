combine_dfs<-function(flist,ftype,rfn="",textfn="",csvfn=""){
  df_data<-data.frame(NULL)
  for (f in flist){
    if (ftype=="R"){
      open(f)
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
    save(list=c("df_data"),file=rfn) 
  }
  if (textfn!=""){
    write.table(df_data,file=textfn,sep="\t",row.names=FALSE,quote=FALSE)
  }
  if (csvfn!=""){
    write.csv(df_tot,file=csvfn,row.names=FALSE,quote=FALSE)
  }
  return(df_data)
}