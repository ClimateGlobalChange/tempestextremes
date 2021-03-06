suppressPackageStartupMessages(require("argparse"))
#Feature requests for BlobStats:
#Line #1: start date/time (as reference for datetime)
#Line #2: name of columns

parser<-ArgumentParser()
#The namelist file provides all of the arguments
#Either use the master namelist file or the individual namelist files
parser$add_argument("-nl","--namelist",help="Master namelist for all flags.",default="")
parser$add_argument("-gn","--generatenl",help="Generate the master namelist",action="store_true")
parser$add_argument("-sl","--setuplist",help="Setup file name file name (mandatory) for --generatenl",default="")
#Most basic function: read BlobStats output to a data table------
parser$add_argument("-rf","--readfiles",help="Tell program to read in BlobStats data",action="store_true")
parser$add_argument("-nlrf","--namelistrf",help="Namelist file name (mandatory) for --readfiles",default="")
#parser$add_argument("-th","--thourly",help="Time resolution (i.e. 6 for 6 hourly)",default="")
#Read in existing file(s) (need to specify type using one of above flags)--------
parser$add_argument("-rt","--readtable",help="Read in existing data file(s) from previous --readfiles output",action="store_true")
parser$add_argument("-nlrt","--namelistrt",help="Namelist file name for --readtable",default="")
#Merge StitchBlobs and DetectBlobs outputs--------
parser$add_argument("-mt","--mergetable",help="Combine outputs from StitchBlobs and DetectBlobs into a table with individual blob information",action="store_true")
parser$add_argument("-nlmt","--namelistmt",help="Namelist file name for --mergetable",default="")
#Make a summary table----------
parser$add_argument("-st","--summarize",help="Make a summary table of the BlobStats data",action="store_true")
parser$add_argument("-nlst","--namelistst",help="Namelist file name for --summarize",default="")
#Read in NetCDFs----------
parser$add_argument("-rn","--readnetcdf",help="Read in NetCDF files",action="store_true")
parser$add_argument("-nlrn","--namelistrn",help="Namelist file name for --readnetcdf",default="")
#Intercomparison of datasets---------
parser$add_argument("-ic","--intercomparison",help="Pearson correlation, probability of co-occurrence between two datasets, and spatial similarity between individual blobs.",action="store_true")
parser$add_argument("-nlic","--namelistic",help="Namelist file name for --intercomparison",default="")
#Pearson and RMSE----------
parser$add_argument("-pr","--pearsonrmse",help="Pearson correlation and RMSE",action="store_true")
parser$add_argument("-nlpr","--namelistpr",help="Namelist file name for --pearsonrmse",default="")
#Generate report
parser$add_argument("-gr","--genreport",help="Generate a summary report on blocking data for specified datasets",action="store_true")
parser$add_argument("-nlgr","--namelistgr",help="namelist file name for --genreport",default="")
#NOW PARSE ALL THE ARGUMENTS--------------
args<-parser$parse_args()

parse_namelist<-function(var,i){
  var_out<-ifelse(length(var)==1,var,var[i])
  return(var_out)
}

##Generate the namelist for running all of the calculations-------------
if (args$generatenl){
  if (args$setuplist==""){
    stop("Provide the name of the setup list file.")
  }
  
  source(args$setuplist)
  setwd(work_dir)
  print("sourcing GenerateNamelist")
  source("generateNamelist.R")
}

#MOST BASIC OPTION: PARSE BLOBSTATS---------
#This will read input from StitchBlobs, input from DetectBlobs, or both
if (args$readfiles){
  if (args$namelistrf=="" & args$namelist==""){
    stop("Must provide the namelist file for --readfiles.")
  }
  
  nl_readfiles<-ifelse(args$namelistrf!="",args$namelistrf,args$namelist)
  source(nl_readfiles)
  setwd(work_dir)
  
  source("readfiles.R")

  for (i in 1:nrun_rf){
    

    varname_i<-parse_namelist(varname,i)
    filename_stitchblobs_i<-parse_namelist(filename_stitchblobs,i)
    filelist_stitchblobs_i<-parse_namelist(filelist_stitchblobs,i)
    rfn_stitch_i<-parse_namelist(rfn_stitch,i)
    df_stitchname_i<-parse_namelist(df_stitchname,i)
    txt_stitch_i<-parse_namelist(txt_stitch,i)
    csv_stitch_i<-parse_namelist(csv_stitch,i)
    
    if (filename_stitchblobs_i=="" & filelist_stitchblobs_i==""){
      stop("Need to either provide a filename (filename_stitchblobs) or list of filenames (filelist_stitchblobs) in the namelist.")
    }
    
    #Vector to hold the BlobStats filenames
    blobstats_vec<-c()
    #Single file
    if (filename_stitchblobs_i!=""){
      blobstats_vec<-c(filename_stitchblobs_i)
    }
    #List of files
    if (filelist_stitchblobs_i!=""){
      fns<-readLines(filelist_stitchblobs_i)
      blobstats_vec<-c(blobstats_vec,fns)
    }
    df_stitch<-read_stats_to_table(blobstats_vec,varname_i,
                                   rfn_stitch_i,txt_stitch_i,csv_stitch_i,
                                   df_stitchname_i)
  }

  
}


#READ IN EXISTING DATA FROM PREVIOUS SESSION (FORMATTED USING ABOVE TOOL)-----
if (args$readtable){
  
  if (args$namelistrt=="" & args$namelist==""){
    stop("Must provide the namelist file for --readtable.")
  }
  nl_readtable<-ifelse(args$namelistrt!="",args$namelistrt,args$namelist)
  source(nl_readtable)
  setwd(work_dir)
  source("readtable.R")
  
  for (i in 1:nrun_rt){
    ftype_rt_i=parse_namelist(ftype_rt,i)
    filename_read_i<-parse_namelist(filename_read,i)
    filelist_read_i<-parse_namelist(filelist_read,i)
    rfn_combine_i<-parse_namelist(rfn_combine,i)
    df_combinename_i<-parse_namelist(df_combinename,i)
    txt_combine_i<-parse_namelist(txt_combine,i)
    csv_combine_i<-parse_namelist(csv_combine,i)
    
    if (filename_read_i=="" & filelist_read_i==""){
      stop("Need to either provide a filename (filename_read) or list of filenames (filelist_read) in the namelist.")
    }
    
    dat_vec<-c()
    #Single file
    if (filename_read_i!=""){
      dat_vec<-c(filename_read_i)
    }
    #List of files
    if (filelist_read_i!=""){
      fns<-readLines(filelist_read_i)
      dat_vec<-c(dat_vec,fns)
    }
    
    
    #Read all files into a single large dataframe
    df_data<-combine_dfs(dat_vec,ftype_rt_i,
                         rfn_combine_i,txt_combine_i,csv_combine_i,
                         df_combinename_i)
  }
  
}

#MERGE STITCHBLOBS AND DETECTBLOBS OUTPUTS----------
if (args$mergetable){
  if (args$namelistmt=="" & args$namelist==""){
    stop("Must provide the namelist file for --mergetable.")
  }
  
  nl_mergetable<-ifelse(args$namelistmt!="",args$namelistmt,args$namelist)
  source(nl_mergetable)
  setwd(work_dir)
  source("readtable.R")
  source("mergetable.R")
  for (i in 1:nrun_mt){
    ftype_mt_i=parse_namelist(ftype_mt,i)
    stitch_file_i<-parse_namelist(stitch_file,i)
    detect_file_i<-parse_namelist(detect_file,i)
    stitch_list_i<-parse_namelist(stitch_list,i)
    detect_list_i<-parse_namelist(detect_list,i)
    rfn_merged_i<-parse_namelist(rfn_merged,i)
    df_merged_name_i<-parse_namelist(df_merged_name,i)
    txt_merged_i<-parse_namelist(txt_merged,i)
    csv_merged_i<-parse_namelist(csv_merged,i)
    if (stitch_file_i=="" & stitch_list_i=="" ){
      stop("Need to provide input StitchBlobs data.")
    }
    if (detect_file_i=="" & detect_list_i==""){
      stop("Need to provide input DetectBlobs data.")
    }
    
    #Make a list of the StitchBlobs files
    stitchvec<-c()
    if (stitch_file_i!=""){
      stitchvec<-c(stitch_file_i)
    }
    if (stitch_list_i!=""){
      fns<-readLines(stitch_list_i)
      stitchvec<-c(stitchvec,fns)
    }
    detectvec<-c()
    if (detect_file_i!=""){
      detectvec<-c(detect_file_i)
    }
    if (detect_list_i!=""){
      fns<-readLines(detect_list_i)
      detectvec<-c(detectvec,fns)
    }
    if (length(stitchvec)<1 | length(detectvec)<1){
      stop("Need to provide files from both StitchBlobs and DetectBlobs.")
    }
    
    #Use the readtable function to load the data
    df_stitch<-combine_dfs(stitchvec,ftype_mt_i)
    df_detect<-combine_dfs(detectvec,ftype_mt_i)
    
    df_merged<-merge_dfs(df_stitch,df_detect,rfn_merged_i,txt_merged_i,csv_merged_i,df_merged_name_i)
    
  }
}

#DATA PROCESSING: GENERATE SUMMARY TABLE-----
if (args$summarize){
  if (args$namelistst=="" & args$namelist==""){
    stop("Must provide the namelist file for --summarize.")
  }
  
  nl_summarize<-ifelse(args$namelistst!="",args$namelistst,args$namelist)
  source(nl_summarize)
  setwd(work_dir)
  source("readtable.R")
  source("summarize.R")
  
  for (i in 1:nrun_st){
    df_input_summ<-data.frame(NULL)
    
    nrun_st_i<-parse_namelist(nrun_st,i)
    ftype_st_i<-parse_namelist(ftype_st,i)
    filename_summ_i<-parse_namelist(filename_summ,i)
    filelist_summ_i<-parse_namelist(filelist_summ,i)
    keepmerge_i<-parse_namelist(keepmerge,i)
    rfn_summ_i<-parse_namelist(rfn_summ,i)
    df_summ_name_i<-parse_namelist(df_summ_name,i)
    txt_summ_i<-parse_namelist(txt_summ,i)
    csv_summ_i<-parse_namelist(csv_summ,i)
    

      #Read in the existing data using combine_tables
      #Create the file list
      if (filename_summ_i=="" & filelist_summ_i==""){
        stop("Need to provide either a filename (filename_summ) or a file list (filelist_summ) for --summarize.")
      }
      dat_vec<-c()
      #Single file
      if (filename_summ_i!=""){
        dat_vec<-c(filename_summ_i)
      }
      #List of files
      if (filelist_summ_i!=""){
        fns<-readLines(filelist_summ_i)
        dat_vec<-c(dat_vec,fns)
      }
      

      df_input_summ<-combine_dfs(dat_vec,ftype_st_i)
      df_summ<-gen_summary_table(df_input_summ,keepmerge_i,
                                 rfn_summ_i,txt_summ_i,csv_summ_i,df_summ_name_i)
  }

}

#READ IN NETCDF DATA----------
if (args$readnetcdf){
  if (args$namelistrn=="" & args$namelist==""){
    stop("Must provide the namelist file for --readnetcdf.")
  }
  
  nl_readnetcdf<-ifelse(args$namelistrn!="",args$namelistrn,args$namelist)
  source(nl_readnetcdf)
  setwd(work_dir)
  source("readnetcdf.R")
  for (i in 1:nrun_rn){
    filename_netcdf_i<-parse_namelist(filename_netcdf,i)
    filelist_netcdf_i<-parse_namelist(filelist_netcdf,i)
    varvec_i<-unlist(parse_namelist(varvec,i))
    outvec_i<-unlist(parse_namelist(outvec,i))
    outrdata_i<-parse_namelist(outrdata,i)
    outnetcdf_i<-parse_namelist(outnetcdf,i)
    timename_i<-parse_namelist(timename,i)
    levname_i<-parse_namelist(levname,i)
    latname_i<-parse_namelist(latname,i)
    lonname_i<-parse_namelist(lonname,i)
    minlat_i<-parse_namelist(minlat,i)
    maxlat_i<-parse_namelist(maxlat,i)
    minlon_i<-parse_namelist(minlon,i)
    maxlon_i<-parse_namelist(maxlon,i)
    minlev_i<-parse_namelist(minlev,i)
    maxlev_i<-parse_namelist(maxlev,i)
    transformto180_i<-parse_namelist(transformto180,i)
    transformto360_i<-parse_namelist(transformto360,i)
    regridto1degree_i<-parse_namelist(regridto1degree,i)

    if (filename_netcdf_i=="" & filelist_netcdf_i==""){
      stop("Need to either provide a filename (filename_netcdf) or list of filenames (filelist_netcdf).")
    }
    
    dat_vec<-c()
    #Single file
    if (filename_netcdf_i!=""){
      dat_vec<-c(filename_netcdf_i)
    }
    #List of files
    if (filelist_netcdf_i!=""){
      print(sprintf("Reading list %s",filelist_netcdf_i))
      fns<-readLines(filelist_netcdf_i)
      dat_vec<-c(dat_vec,fns)
    }
    
    if (length(varvec_i)!=length(outvec_i)){
      stop("Length of variable lists for inputs and outputs differ. Check input and output variable list strings.")
    }
    
    vars_list<-read_netcdf(dat_vec,varvec_i,outvec_i,timename_i,levname_i,latname_i,lonname_i,
                           minlat_i,maxlat_i,minlon_i,maxlon_i,
                           minlev_i,maxlev_i,outnetcdf_i,outrdata_i,
                           transformto180_i,transformto360_i,regridto1degree_i)
    
  }
  
}

#COMPARE TWO DATASETS--------------------------------
if (args$intercomparison){
  if (args$namelistic=="" & args$namelist==""){
    stop("Must provide the namelist file for --intercomparison.")
  }
  
  nl_prob<-ifelse(args$namelistic!="",args$namelistic,args$namelist)
  source(nl_prob)
  setwd(work_dir)
  source("intercomparison.R")
  for (i in 1:nrun_ic){
    table_file_1_i<-parse_namelist(table_file_1,i)
    table_file_2_i<-parse_namelist(table_file_2,i)
    df_name_1_i<-parse_namelist(df_name_1,i)
    df_name_2_i<-parse_namelist(df_name_2,i)
    blob_file_1_i<-parse_namelist(blob_file_1,i)
    blob_file_2_i<-parse_namelist(blob_file_2,i)
    var_name_1_i<-parse_namelist(var_name_1,i)
    var_name_2_i<-parse_namelist(var_name_2,i)
    regrid_i<-parse_namelist(regrid,i)
    rfn_ps_i<-parse_namelist(rfn_ps,i)
    txt_overlaps_i<-parse_namelist(txt_overlaps,i)
    txt_ps_i<-parse_namelist(txt_ps,i)
    #Load the files-- careful not to overwrite!
    load(table_file_1_i)
    df1<-get(df_name_1_i)
    df2<-NULL
    if (table_file_2_i!=""){
      load(table_file_2_i)
      df2<-get(df_name_2_i)
    }
    #Load the blob files
    #File 1
    load(blob_file_1_i)
    blob1<-get(var_name_1_i)
    lat1<-lat_axis
    lon1<-lon_axis
    time1<-time_format
    #File 2
    load(blob_file_2_i)
    blob2<-get(var_name_2_i)
    lat2<-lat_axis
    lon2<-lon_axis
    time2<-time_format
    
    probsim<-overlaps_calc(df1,blob1,time1,lat1,lon1,blob2,
                           df2,time2,lat2,lon2,rfn_ps_i,txt_overlaps_i,txt_ps_i,regrid_i)
    
  }
  
}

if (args$pearsonrmse){
  if (args$namelistpr=="" & args$namelist==""){
    stop("Must provide the namelist file for --pearsonrmse.")
  }
  
  nl_prob<-ifelse(args$namelistpr!="",args$namelistpr,args$namelist)
  source(nl_prob)
  setwd(work_dir)
  source("pearsonrmse.R")
  for (i in 1:nrun_pr){
    blob_file_1_i<-parse_namelist(blob_file_1,i)
    blob_file_2_i<-parse_namelist(blob_file_2,i)
    data_name_1_i<-parse_namelist(data_name_1,i)
    data_name_2_i<-parse_namelist(data_name_2,i)
    var_name_1_i<-parse_namelist(var_name_1,i)
    var_name_2_i<-parse_namelist(var_name_2,i)
    regrid_i<-parse_namelist(regrid,i)
    rfn_pr_i<-parse_namelist(rfn_pr,i)
    txt_pr_i<-parse_namelist(txt_pr,i)
    #Load the blob files
    #File 1
    load(blob_file_1_i)
    blob1<-get(var_name_1_i)
    lat1<-lat_axis
    lon1<-lon_axis
    time1<-time_format
    #File 2
    load(blob_file_2_i)
    blob2<-get(var_name_2_i)
    lat2<-lat_axis
    lon2<-lon_axis
    time2<-time_format
    
    pr<-pearson_rmse(data_name_1_i,data_name_2_i,blob1,time1,lat1,lon1,
                     blob2,time2,lat2,lon2,regrid_i,useCommonTime,
                     rfn_pr_i,txt_pr_i)
    
  }
  
}



#GENERATE A REPORT WITH ALL OF THE PROVIDED DATA
if (args$genreport){
  if (args$namelistgr=="" & args$namelist==""){
    stop("Must provide the namelist file for --genreport.")
  }
  nl_report<-ifelse(args$namelistgr!="",args$namelistgr,args$namelist)
  source(nl_report)
  setwd(work_dir)
  source("generateReport.R")
}







