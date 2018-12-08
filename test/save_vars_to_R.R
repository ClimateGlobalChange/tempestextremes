#!/usr/bin/env Rscript

#This code will read a variable to an RData file
#This is being written in order to clear out some of the NetCDF files
# that are taking up a lot of space in Agri!

library(RNetCDF)
library(argparse)

source("~/tempestextremes/test/sector_funcs.R")

parser<- ArgumentParser()

#Create the flags
parser$add_argument("-fn","--filename",help="Name of input netCDF file",default="")
parser$add_argument("-fl","--filelist",help="Name of input netCDF files list",default="")
parser$add_argument("-vn","--varname",help="Variable name to be stored",default="")
parser$add_argument("-vl","--varlist",help="String of variable names to be stored (separated by commas, no space, example: T,U,V)",default="")
parser$add_argument("-ovl","--outvarlist",help="String of output array names (separated by commas, no space, example: T,U,V",default="")
parser$add_argument("-o","--out",help="Name of output RData file",default="")
parser$add_argument("--region",help="Subset region for input data",default="")
parser$add_argument("--syear",default=0,help="Starting year")
parser$add_argument("--eyear",default=0,help="Ending year")
parser$add_argument("--rev_lat",help="Flip the orientation of the latitude axis",default=FALSE)
parser$add_argument("--lon360",help="Longitude axis goes from 0->360 rather than -180->180",default=FALSE)
parser$add_argument("--latname",default="lat")
parser$add_argument("--lonname",default="lon")
parser$add_argument("--nperday",default=4)

args<-parser$parse_args()
if (args$filename=="" & args$filelist==""){
  stop("Need to either provide a filename or list of filenames.")
}
if (args$varname=="" & args$varlist==""){
  stop("Need to provide at least one variable name.")
}
if (args$region==""){
  stop("Need to specify the region (NA, NC, NP, SA, SI, SP)")
}
if (args$syear==0 | args$eyear==0){
  stop("Need to provide start and end years.")
}

#The regions and longitude ranges
regions<-c("NA","NC","NP","SA","SI","SP")
MIN_LAT=c(25, 25, 25, -75, -75, -75)
MAX_LAT=c(75, 75, 75, -25, -25, -25)
if (args$lon360==TRUE){
  LEFT_BOUND=c(250, 30, 130, 290, 20, 120)
  RIGHT_BOUND=c(50, 150, 270, 40, 140, 310)
}else{
  LEFT_BOUND=c(-110, 30, 130, -70, 20, 120)
  RIGHT_BOUND=c(50, 150, -90, 40, 140, -50)
}

#Create the vector of files
fname_list<-c()
if (args$filename!=""){
  fname_list[1]<-args$filename
}else{
  #Parse the list text file
  fname_list<-readLines(args$filelist)
}

#Create the vector of variable names
vname_list<-c()
if (args$varname!=""){
  vname_list[1]<-args$varname
}else{
  vname_list<-unlist(strsplit(args$varlist,split=","))
}
vstr=paste(vname_list,collapse=", ")
print(sprintf("Reading in variables %s",vstr))

if (args$outvarlist!=""){
  outvlist<-unlist(strsplit(args$outvarlist,split=","))
}else{
  outvlist<-vname_list
}

#Open the first file and get the axes for lat, lon
nf<-open.nc(fname_list[1])
lat<-var.get.nc(nf,args$latname)
nlat<-length(lat)
lon<-var.get.nc(nf,args$lonname)
nlon<-length(nlon)
close.nc(nf)

#Get the region specific info
seci<-which(regions==args$region)
lb<-get_index(lon,LEFT_BOUND[seci])
rb<-get_index(lon,RIGHT_BOUND[seci])
if (lb>rb){
  lons_sub<-c(seq(lb,length(lons)),seq(1,rb))
}else{
  lons_sub<-seq(lb,rb)
}
tb<-get_index(lat,MIN_LAT[seci])
bb<-get_index(lat,MAX_LAT[seci])
if (tb>bb){
  lats_sub<-seq(bb,tb)
}else{
  lats_sub<-seq(tb,bb)
}
latsize<-length(lats_sub)
lonsize<-length(lons_sub)
nyears=args$eyear-args$syear+1

lats_seq<-lat(lats_sub)
lons_seq<-lon(lons_sub)

tsz<-args$nperday*nyears*365

tempFrame<-array(NA,c(lonsize,latsize,tsz))
#Create the dataframes
for (v in vname_list){
  
}


for (f in fname_list){
  nf<-open.nc(f){
    #Read in the latitude variable
    #Now read each variable in
  }
  close.nc(f)
}