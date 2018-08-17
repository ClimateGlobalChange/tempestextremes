suppressPackageStartupMessages(require("argparse"))
#Feature requests for BlobStats:
#Line #1: start date/time (as reference for datetime)
#Line #2: name of columns

parser<-ArgumentParser()
#Takes either a list or a single file
#parser$add_argument("-v","--verbose",help="Print each table to the screen.",action="store_true")
parser$add_argument("-fn","--filename",help="Name of input file",default="")
parser$add_argument("-fl","--filelist",help="Name of input list of files",default="")
#Most basic function: read BlobStats output to a data table------
parser$add_argument("-rf","--readfiles",help="Tell program to read in BlobStats data",action="store_true")
#Currently: manually specify time resolution
parser$add_argument("-th","--thourly",help="Time resolution (i.e. 6 for 6 hourly)",default="")
parser$add_argument("-an","--algname",help="Type of blocking algorithm used",default="")
parser$add_argument("--rtable",help="Name of output RData file for --readfiles or --readtable",default="")
parser$add_argument("--texttable",help="Name of output text file for --readfiles or --readtable",default="")
parser$add_argument("--csvtable",help="Name of output CSV file for --readfiles or --readtable",default="")

#Read in existing file(s) (need to specify type using one of above flags)--------
parser$add_argument("-rt","--readtable",help="Read in existing data file(s) from previous --readfiles output",action="store_true")
parser$add_argument("--isr",help="Tells --readtable that input is in RData format",action="store_true")
parser$add_argument("--istext",help="Tells --readtable that input is in text format",action="store_true")
parser$add_argument("--iscsv",help="Tells --readtable that input is in CSV format",action="store_true")
#Make a summary table----------
parser$add_argument("-st","--summarize",help="Make a summary table of the BlobStats data",action="store_true")
parser$add_argument("--rsumm",help="Name of output RData file for --summarize",default="")
parser$add_argument("--textsumm",help="Name of output text file for --summarize",default="")
parser$add_argument("--csvsumm",help="Name of output CSV file for --summarize",default="")
#Read in NetCDFs----------
parser$add_argument("-rn","--readnetcdf",help="Read in NetCDF files",action="store_true")
parser$add_argument("-vl","--varlist",help="List of input NetCDF variable(s) in format VAR1,VAR2,VAR3 (no spaces between commas)",default="")
parser$add_argument("-ovl","--outvarlist",help="Optional list of output variable names corresponding to --varlist inputs",default="")
parser$add_argument("--outnetcdf",help="Name of output NetCDF file for --readnetcdf",default="")
parser$add_argument("--outrdata",help="Name of output RData file for --readnetcdf",default="")
parser$add_argument("--subset",help="Subset the lat/lon extent",action="store_true")
parser$add_argument("--minlat",help="Minimum latitude value",default=NULL)
parser$add_argument("--maxlat",help="Maximum latitude value",default=NULL)
parser$add_argument("--minlon",help="Leftmost longitude value",default=NULL)
parser$add_argument("--maxlon",help="Rightmost longitude value",default=NULL)
parser$add_argument("--minlev",help="Minimum level value",default=NULL)
parser$add_argument("--maxlev",help="Maximum level value",default=NULL)
parser$add_argument("--timename",help="Name of NetCDF time variable",default="time")
parser$add_argument("--levname",help="Name of NetCDF level variable",default="lev")
parser$add_argument("--latname",help="Name of NetCDF latitude variable",default="lat")
parser$add_argument("--lonname",help="Name of NetCDF longitude variable",default="lon")
#Combine Rdata files into a single Rdata file
#parser$add_argument("-cd","--combinedata",help="Combine RData files into a single output RData file",action="store_true")

#Intercomparison of datasets---------
parser$add_argument("-pr","--probability",help="Probability of co-occurrence between two .",default="")

#NOW PARSE ALL THE ARGUMENTS--------------
args<-parser$parse_args()

#Data frames for holding the info
df_data<-data.frame(NULL)
df_summ<-data.frame(NULL)
list_vars<-list()
#MOST BASIC OPTION: PARSE BLOBSTATS---------
if (args$readfiles){
  source("read_stitch.R")
  if (args$filename=="" & args$filelist==""){
    stop("Need to either provide a filename (-fn) or list of filenames (-fl).")
  }
  if (args$thourly==""){
    stop("Need to provide time resolution (-th).")
  }

  #Vector to hold the BlobStats filenames
  blobstats_vec<-c()
  #Number of time steps
  nhrs<-as.numeric(args$thourly)

  #Single file
  if (args$filename!=""){
    blobstats_vec<-c(args$filename)
  }
  #List of files
  if (args$filelist!=""){
    fns<-readLines(args$filelist)
    blobstats_vec<-c(blobstats_vec,fns)
  }
  df_data<-read_stats_to_table(blobstats_vec,nhrs,args$algname,
                               args$rtable,args$texttable,args$csvtable)
  if (args$verbose){
    print(df_data)
  }
}

#READ IN EXISTING DATA FROM PREVIOUS SESSION (FORMATTED USING ABOVE TOOL)-----
if (args$readtable){
  if (!args$rfile && !args$textfile && !args$csvfile){
    stop("Need to specify the format of the input data: RData (--rfile), text (--textfile), or CSV (--csvfile)")
  }
  if (args$filename=="" & args$filelist==""){
    stop("Need to either provide a filename (-fn) or list of filenames (-fl).")
  }
  dat_vec<-c()
  #Single file
  if (args$filename!=""){
    dat_vec<-c(args$filename)
  }
  #List of files
  if (args$filelist!=""){
    fns<-readLines(args$filelist)
    dat_vec<-c(dat_vec,fns)
  }
  
  ftype=""
  if (args$isr){
    ftype="R"
  }else if (args$iscsv){
    ftype="CSV"
  }else if (args$istext){
    ftype="text"
  }
  #Read all files into a single large dataframe
  source("combine_tables.R")
  df_data<-combine_dfs(dat_vec,ftype,args$rtable,args$texttable,args$csvtable)

}
#print(ncol(df_data))
#DATA PROCESSING: GENERATE SUMMARY TABLE-----
if (args$summarize){
  source("summ_table.R")
  if (ncol(df_data)<1){
    stop("Need to either read in BlobStats data (-rf) or an existing table (-rt)")
  }
  df_summ<-gen_summary_table(df_data,args$rsumm,args$textsumm,args$csvsumm)
}
#READ IN NETCDFS FOR BLOBS, Z500, ETC-------
#Either saves to RData or NetCDF
if (args$readnetcdf){
  if (args$filename=="" & args$filelist==""){
    stop("Need to specify either a file name (-fn) or file list (-fl) of netCDF files to read in.")
  }

  nlist<-c()
  if (args$filename!=""){
    nlist<-c(nlist,args$filename)
  }
  if (args$filelist!=""){
    fns<-readLines(args$filelist)
    nlist<-c(nlist,fns)
  }
  if (args$varlist==""){
    stop("Need to specify the variable names of interest (-vl).")
  }
  vvec<-unlist(strsplit(args$varlist,split=","))
  if (args$outvarlist!=""){
    ovec<-unlist(strsplit(args$outvarlist,split=","))
    if (length(vvec)!=length(ovec)){
      stop("The length of the two variables lists does not match. Check --varlist and --outvarlist inputs.")
    }
  }else{
    ovec<-vvcec
  }
  #Read the NetCDF files to the workspace (or save to RData)
  source("read_netCDF_to_R.R")
  list_vars<-read_netcdf(nlist,vvec,ovec,args$timename,args$levname,args$latname,args$lonname,args$subset,args$minlat,args$maxlat,
              args$minlon,args$maxlon,args$minlev,args$maxlev,args$outnetcdf,args$outrdata)
}


