#This will generate a master namelist based on the inputs 
# from the report namelist
dir_check<-sprintf("%s/%s",output_dir,output_subdir)
ifelse(!dir.exists(file.path(dir_check)),dir.create(file.path(dir_check)),FALSE)
fname_namelist_out<-sprintf("%s/%s",dir_check,fname_namelist)
#READFILES
if (!use_detectblob){
  flist<-sprintf("%s",stitch_lists)
  var_inputs<-Varnames
  suffix_table<-rep("stitchtable",length(stitch_lists))
}else{
  flist<-sprintf("%s",c(stitch_lists,detect_lists))
  var_inputs<-c(Varnames,Varnames)
  suff1<-rep("stitchtable",length(stitch_lists))
  suff2<-rep("nostitchtable",length(detect_lists))
  suffix_table<-c(suff1,suff2)
}
nrun_rf<-length(flist)
#Will always output RData files
rdata_readname<-sprintf("%s/%s/%s_%s_%s.RData",output_dir,output_subdir,
                        output_prefix,var_inputs,suffix_table)
#optional txt/csv files
txt_readname<-ifelse(output_csv==TRUE,sprintf("%s/%s/%s_%s_%s.txt",
                                               output_dir,output_subdir,
                                               output_prefix,var_inputs,
                                               suffix_table),"")
csv_readname<-ifelse(output_csv==TRUE,sprintf("%s/%s/%s_%s_%s.csv",
                                               output_dir,output_subdir,
                                               output_prefix,var_inputs,
                                               suffix_table),"")
#collapsed strings for output
var_string<-paste(var_inputs,collapse="\",\"")
res_string<-ifelse(!exists("resolutions"),"",paste(resolutions,collapse="\",\""))
flist_string<-paste(flist,collapse="\",\"")
rdata_readstring<-paste(rdata_readname,collapse="\",\"")
csv_readstring<-ifelse(csv_readname=="","\"\"",paste(csv_readname,collapse="\",\""))
txt_readstring<-ifelse(txt_readname=="","\"\"",paste(txt_readname,collapse="\",\""))
if (use_detectblob){
  detect_rfiles<-rdata_readname[(length(stitch_lists)+1):length(flist)]
  detect_string<-paste(detect_rfiles,collapse="\",\"")
}
stitch_rfiles<-rdata_readname[1:length(stitch_lists)]
stitch_string<-paste(stitch_rfiles,collapse="\",\"")
#MERGETABLE
if (use_detectblob){
  #Filenames for mergetable
  rdata_mergenames<-sprintf("%s/%s/%s_%s_mergetable.RData",output_dir,output_subdir,
                            output_prefix,Varnames)
  rdata_mergestring<-paste(rdata_mergenames,collapse="\",\"")
  
  txt_mergename<-ifelse(output_csv==TRUE,sprintf("%s/%s/%s_%s_mergetable.txt",
                                                output_dir,output_subdir,
                                                output_prefix,Varnames),"")
  csv_mergename<-ifelse(output_csv==TRUE,sprintf("%s/%s/%s_%s_mergetable.csv",
                                                output_dir,output_subdir,
                                                output_prefix,Varnames),"")
  csv_mergestring<-ifelse(csv_mergename=="","\"\"",paste(csv_mergename,collapse="\",\""))
  txt_mergestring<-ifelse(txt_mergename=="","\"\"",paste(txt_mergename,collapse="\",\""))
}
#SUMMTABLE
#USES BLOBSTATS TABLE IF MERGETABLE IS NOT RUN
if (use_detectblob){
  summ_inputs<-rdata_mergenames
}else{
  summ_inputs<-stitch_rfiles
}
summ_string<-paste(summ_inputs,collapse="\",\"")
rdata_summ<-sprintf("%s/%s/%s_%s_summtable.RData",output_dir,output_subdir,
                    output_prefix,Varnames)
rdata_summstring<-paste(rdata_summ,collapse="\",\"")

txt_summname<-ifelse(output_csv==TRUE,sprintf("%s/%s/%s_%s_summtable.txt",
                                               output_dir,output_subdir,
                                               output_prefix,Varnames),"")
csv_summname<-ifelse(output_csv==TRUE,sprintf("%s/%s/%s_%s_summtable.csv",
                                               output_dir,output_subdir,
                                               output_prefix,Varnames),"")
csv_summstring<-ifelse(csv_summname=="","\"\"",paste(csv_summname,collapse="\",\""))
txt_summstring<-ifelse(txt_summname=="","\"\"",paste(txt_summname,collapse="\",\""))

#READNETCDF
#Collapse the netcdf list 
stitchblob_string<-paste(sprintf("%s",stitchblob_lists),collapse="\",\"")
varvec_string<-paste(varvec,collapse="\",\"")
rdata_blobname<-sprintf("%s/%s/%s_%s_blobdata.RData",output_dir,output_subdir,
                        output_prefix,Varnames)
rdata_blobstring<-paste(rdata_blobname,collapse="\",\"")
TF180<-ifelse(transformto180==TRUE,"TRUE","FALSE")
TF360<-ifelse(transformto360==TRUE,"TRUE","FALSE")
TFregrid<-ifelse(regridto1degree==TRUE,"TRUE","FALSE")
#INTERCOMPARISON
nvars<-length(Varnames)

table1<-c()
table2<-c()
blob1<-c()
blob2<-c()
var1<-c()
var2<-c()
icnames<-c()
txtico<-c()
txticp<-c()

suff<-ifelse(use_detectblob==TRUE,"mergetable.RData","stitchtable.RData")
for (i in 1:length(Varnames)){
  for (j in 2:length(Varnames)){
    if (i<j){
      t1<-sprintf("%s/%s/%s_%s_%s",output_dir,output_subdir,
                  output_prefix,Varnames[i],suff)
      t2<-sprintf("%s/%s/%s_%s_%s",output_dir,output_subdir,
                  output_prefix,Varnames[j],suff)
      b1<-sprintf("%s/%s/%s_%s_blobdata.RData",output_dir,output_subdir,
                  output_prefix,Varnames[i])
      b2<-sprintf("%s/%s/%s_%s_blobdata.RData",output_dir,output_subdir,
                  output_prefix,Varnames[j])
      iname<-sprintf("%s/%s/%s_ic_%s_%s.RData",output_dir,output_subdir,
                     output_prefix,Varnames[i],Varnames[j])
      iconame<-sprintf("%s/%s/%s_ic_%s_%s_table.txt",output_dir,output_subdir,
                       output_prefix,Varnames[i],Varnames[j])
      icpname<-sprintf("%s/%s/%s_ic_%s_%s_probsim.txt",output_dir,output_subdir,
                       output_prefix,Varnames[i],Varnames[j])
      table1<-c(table1,t1)
      table2<-c(table2,t2)
      blob1<-c(blob1,b1)
      blob2<-c(blob2,b2)
      var1<-c(var1,varvec[i])
      var2<-c(var2,varvec[j])
      icnames<-c(icnames,iname)
      txtico<-c(txtico,iconame)
      txticp<-c(txticp,icpname)
    }
  }
}

numcombos<-length(blob1)
table1_string<-paste(table1,collapse="\",\n\"")
table2_string<-paste(table2,collapse="\",\n\"")
blob1_string<-paste(blob1,collapse="\",\n\"")
blob2_string<-paste(blob2,collapse="\",\n\"")
var1_string<-paste(var1,collapse="\",\n\"")
var2_string<-paste(var2,collapse="\",\n\"")
icname_string<-paste(icnames,collapse="\",\n\"")
txtico_string<-paste(txtico,collapse="\",\n\"")
txticp_string<-paste(txticp,collapse="\",\n\"")
df_icname<-ifelse(use_detectblob==TRUE,"df_merged","df_blob_per_timestep")

#Generating the report!!
Varnames_sep<-paste(Varnames,collapse="_")
reportname_string<-sprintf("%s/%s/%s_%s_report.html",output_dir,output_subdir,
                           output_prefix,Varnames_sep)
Varnames_string<-paste(Varnames,collapse="\",\"")


#WRITING THE NAMELIST!----------

print(sprintf("Writing namelist %s",fname_namelist_out))

sink(fname_namelist_out,split=T)
cat(sprintf("work_dir<-\"%s\" \n",work_dir))
cat(sprintf("output_dir<-\"%s/%s\" \n",output_dir,output_subdir))
if (res_string!=""){
  cat(sprintf("resolutions<-c(\"%s\")\n\n",res_string))  
}
cat("#READFILES SPECS----------\n")
cat(sprintf("nrun_rf<-%d \n",nrun_rf))
cat(sprintf("varname<-c(\"%s\") \n",var_string))
cat(sprintf("filelist_stitchblobs<-c(\"%s\") \n",flist_string))
cat("filename_stitchblobs<-\"\"\n")
cat(sprintf("rfn_stitch<-c(\"%s\") \n",rdata_readstring))
cat("df_stitchname<-\"df_blob_per_timestep\" \n")
if(txt_readname==""){
  cat("txt_stitch<-\"\" \n")
}else{
  cat(sprintf("txt_stitch<-c(\"%s\")\n",txt_readstring))
}
if(csv_readname==""){
  cat("csv_stitch<-\"\" \n")
}else{
  cat(sprintf("csv_stitch<-c(\"%s\")\n",csv_readstring))
}
cat("\n")
if (use_detectblob){
  cat("#MERGETABLE SPECS----------\n")
  cat(sprintf("nrun_mt<-%d \n",nrun_rf/2))
  cat("ftype_mt<-\"R\"\n")
  cat(sprintf("stitch_file<-c(\"%s\")\n",stitch_string))
  cat(sprintf("detect_file<-c(\"%s\")\n",detect_string))
  cat("stitch_list<-\"\"\n")
  cat("detect_list<-\"\"\n")
  cat(sprintf("rfn_merged<-c(\"%s\")\n",rdata_mergestring))
  cat("df_merged_name<-\"df_merged\"\n")
  if(txt_mergename==""){
    cat("txt_merged<-\"\" \n")
  }else{
    cat(sprintf("txt_merged<-c(\"%s\")\n",txt_mergestring))
  }
  if(csv_readname==""){
    cat("csv_merged<-\"\" \n")
  }else{
    cat(sprintf("csv_merged<-c(\"%s\")\n",csv_mergestring))
  }
  cat("\n")
}

cat("#SUMMTABLE SPECS----------\n")
cat(sprintf("nrun_st<-%d\n",length(Varnames)))
cat("ftype_st<-\"R\"\n")
cat(sprintf("filename_summ<-c(\"%s\")\n",summ_string))
cat(sprintf("filelist_summ<-\"\"\n"))
cat("keepmerge<-TRUE\n")
cat(sprintf("rfn_summ<-c(\"%s\")\n",rdata_summstring))
cat("df_summ_name<-\"df_summ\"\n")
if(txt_summname==""){
  cat("txt_summ<-\"\" \n")
}else{
  cat(sprintf("txt_summ<-c(\"%s\")\n",txt_summstring))
}
if(csv_summname==""){
  cat("csv_summ<-\"\" \n")
}else{
  cat(sprintf("csv_summ<-c(\"%s\")\n",csv_summstring))
}
cat("\n")

cat("#READNETCDF SPECS----------\n")
cat(sprintf("nrun_rn<-%d\n",length(Varnames)))
cat("filename_netcdf<-\"\"\n")
cat(sprintf("filelist_netcdf<-c(\"%s\")\n",stitchblob_string))
cat(sprintf("varvec<-c(\"%s\")\n",varvec_string))
cat("outvec<-varvec\n")
cat(sprintf("outrdata<-c(\"%s\")\n",rdata_blobstring))
cat("outnetcdf<-\"\"\n")
cat(sprintf("timename<-\"%s\"\n",timename))
cat("levname<-\"lev\"\n")
cat(sprintf("latname<-\"%s\"\n",latname))
cat(sprintf("lonname<-\"%s\"\n",lonname))
cat(sprintf("transformto180<-%s\n",TF180))
cat(sprintf("transformto360<-%s\n",TF360))
if (exists("minlat")){
  cat(sprintf("minlat<-%f\n",minlat))  
}else{
  cat("minlat<-\"\"\n")
}
if (exists("maxlat")){
  cat(sprintf("maxlat<-%f\n",maxlat))
}else{
  cat("maxlat<-\"\"\n")
}
if (exists("minlon")){
  cat(sprintf("minlon<-%f\n",minlon))  
}else{
  cat("minlon<-\"\"\n")
}
if (exists("maxlon")){
  cat(sprintf("maxlon<-%f\n",maxlon))  
}else{
  cat("maxlon<-\"\"\n")
}
cat("minlev<-\"\"\n")
cat("maxlev<-\"\"\n")
cat(sprintf("regridto1degree<-%s\n",TFregrid))
cat("\n")

cat("#INTERCOMPARE SPECS----------\n")
cat(sprintf("nrun_ic<-%d\n",numcombos))
cat(sprintf("table_file_1<-c(\"%s\")\n",table1_string))
cat(sprintf("table_file_2<-c(\"%s\")\n",table2_string))
cat(sprintf("df_name_1<-\"%s\"\n",df_icname))
cat(sprintf("df_name_2<-\"%s\"\n",df_icname))
cat("regrid<-FALSE\n")
cat(sprintf("blob_file_1<-c(\"%s\")\n",blob1_string))
cat(sprintf("var_name_1<-c(\"%s\")\n",var1_string))
cat(sprintf("blob_file_2<-c(\"%s\")\n",blob2_string))
cat(sprintf("var_name_2<-c(\"%s\")\n",var2_string))
cat(sprintf("rfn_ps<-c(\"%s\")\n",icname_string))
cat(sprintf("txt_overlaps<-c(\"%s\")\n",txtico_string))
cat(sprintf("txt_ps<-c(\"%s\")\n",txticp_string))
cat("\n")

cat("#REPORT SPECS----------\n")
cat(sprintf("output_name<-\"%s\"\n",reportname_string))
cat(sprintf("metadata_datasets<-\"%s\"\n",metadata_datasets))
cat(sprintf("Varnames<-c(\"%s\")\n",Varnames_string))
cat(sprintf("nyears<-%d\n",nyears))
cat(sprintf("mergefiles<-c(\"%s\")\n",summ_string))
cat(sprintf("summfiles<-c(\"%s\")\n",rdata_summstring))
cat(sprintf("blobfiles<-c(\"%s\")\n",rdata_blobstring))
cat(sprintf("blobname<-c(\"%s\")\n",varvec_string))
cat(sprintf("icfiles<-c(\"%s\")\n",icname_string))
sink()

print(sprintf("Wrote namelist file %s",fname_namelist_out))