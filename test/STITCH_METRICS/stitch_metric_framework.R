require("argparse")

#Feature requests for BlobStats:
#Line #1: start date/time (as reference for datetime)
#Line #2: name of columns

parser<-ArgumentParser()
#Most basic function: read BlobStats output to a data table
#Takes either a list or a single file
parser$add_argument("-rf","--readfiles",help="Tell program to read in BlobStats data",action="store_true")
parser$add_argument("-fn","--filename",help="Name of input file",default="")
parser$add_argument("-fl","--filelist",help="Name of input list of files",default="")
parser$add_argument("-vn","--varname",help="Type of blocking algorithm used",default="")
parser$add_argument("--rtable",help="Name of output RData file for --readfiles",default="")
parser$add_argument("--texttable",help="Name of output text file for --readfiles",default="")
parser$add_argument("--csvtable",help="Name of output CSV file for --readfiles",default="")
#Currently: manually specify time resolution
parser$add_argument("-th","--thourly",help="Time resolution (i.e. 6 for 6 hourly)",default="")
#Read in existing file(s) (need to specify type using one of above flags)
parser$add_argument("-rt","--readtable",help="Read in existing data file(s) from previous --readfiles output",action="store_true")
parser$add_argument("--rfile",help="Tells --readtable that input file is in RData format",action="store_true")
parser$add_argument("--textfile",help="Tells --readtable that input file is in text format",action="store_true")
parser$add_argument("--csvfile",help="Tells --readtable that input file is in CSV format",action="store_true")
#Make a summary table
parser$add_argument("-st","--summarize",help="Make a summary table of the BlobStats data",action="store_true")
parser$add_argument("--rsumm",help="Name of output RData file for --summarize",default="")
parser$add_argument("--textsumm",help="Name of output text file for --summarize",default="")
parser$add_argument("--csvsumm",help="Name of output CSV file for --summarize",default="")


args<-parser$parse_args()
df_data<-data.frame(NULL)
df_summ<-data.frame(NULL)
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
  df_data<-read_stats_to_table(blobstats_vec,nhrs,args$varname,
                               args$rtable,args$texttable,args$csvtable)
  
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
  #Read all files into a single large dataframe
  for (f in dat_vec){
    if (args$rfile){
      open(f)
    }
    if (args$csvfile){
      df_tot<-read.csv(f)
    }
    if (args$textfile){
      df_tot<-read.table(f,header=TRUE,sep="\t")
    }
    df_data<-rbind(df_data,df_tot)
  }
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
