import numpy as np
import xarray as xa
import pandas as pd
import argparse


parser=argparse.ArgumentParser(description="Provide a text file list for the linear regression calculation.")
parser.add_argument("-f","--filelist",required=True,action="store",dest="lsname")
parser.add_argument("-o","--out",required=True,action="store")
parser.add_argument("--timename",default="time")
parser.add_argument("--latname",default="lat")
parser.add_argument("--lonname",default="lon")
parser.add_argument("--varname",default="SMOOTHED_Z500")
results=parser.parse_args()

flist= open(results.lsname,"r")
datafiles=flist.read().splitlines()
nyears=len(datafiles)
yrange=np.arange(0.,nyears)

#This is just a pointer to the dataset
dataset=xa.open_mfdataset(datafiles,decode_times=False,concat_dim=results.timename)
#Read in the data
svar=dataset[results.varname][:]
tvar=dataset[results.timename][:]
latvar=dataset[results.latname][:]
lonvar=dataset[results.lonname][:]
nlat=len(latvar)
nlon=len(lonvar)
ndays=dataset.chunks[results.timename][0]
#Initialize the arrays that will store the values
marray=np.zeros((ndays,nlat,nlon))
barray=np.zeros((ndays,nlat,nlon))
print("Initialized arrays. Calculating linear regression coefficients.")
for t in range(0,ndays):
    #Generate the dates for the linear regression
    #nday=int(dataset.time.dt.day[t])
    #nmonth=int(dataset.time.dt.month[t])
    #nyear=int(dataset.time.dt.year[t])
    #dstring="{:04d}/{:02d}/{:02d}".format(nyear,nmonth,nday)
    #print(dstring)
    #tsel=pd.date_range(start=dstring,periods=nyears,freq=pd.DateOffset(years=1))
    #Subset the data by the time steps
    zpts=svar.sel(time=t).values
    print("Subsetting for time {:d}".format(t))
    #Fit the linear regression to each lat/lon
    for x in range(0,nlat):
        for y in range(0,nlon):
            zsub=zpts[:,x,y]
            #Check for nan
            idx=np.isfinite(yrange) & np.isfinite(zsub)
            m,b= np.polyfit(yrange[idx],zsub[idx],1)
            marray[t,x,y]=m
            barray[t,x,y]=b

#Write the output file
print("Writing {} to file".format(results.out))
daysSince=np.arange(0,ndays)
dataset_out=xa.Dataset({'slope':([results.timename,results.latname,results.lonname],marray),
    'intercept':([results.timename,results.latname,results.lonname],barray)},
    coords={results.timename:daysSince,
        results.latname:latvar,
        results.lonname:lonvar})

netcdf_calendar = "standard"
if (ndays<365):
    netcdf_calendar="360_day"

dataset_out[results.timename].attrs={"units":"days since 0001-01-01","calendar":netcdf_calendar}
dataset_out.to_netcdf(results.out)
