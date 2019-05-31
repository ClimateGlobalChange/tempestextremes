import argparse
import numpy as np
import xarray as xa
import pandas as pd
from datetime import datetime
import cftime
import math

def deg2rad(deg):
    return(deg*(math.pi/180.))

def getDistanceFromLatLonInKm(lat1,lon1,lat2,lon2):
    R=6371.
    dlat=deg2rad(lat2-lat1)
    dlon=deg2rad(lon2-lon1)
    a=math.sin(dlat/2.)*math.sin(dlat/2.) + math.cos(deg2rad(lat1)) *\
    math.cos(deg2rad(lat2)) *math.sin(dlon/2.)*math.sin(dlon/2.)
    
    b = 2. * math.atan2(math.sqrt(a),math.sqrt(1-a))
    d=R*b
    return(d)

#Everything is on the 0->360 lon range! 
def calcAI(anom_mask,r_info):
    #print(anom_mask.shape)
    xlen=len(anom_mask.lon)
    ylen=len(anom_mask.lat)
    avals=anom_mask.values
    lons=anom_mask.lon.values
    lats=anom_mask.lat.values
    store_anom=r_info
    anom_vals=[]
    area_vals=[]
    amax = np.max(anom_mask.values)
    for i in range(0,(xlen-1)):
        i1=i
        i2=i+1
        if (ylen<2):
            asub=avals[:,i1:(i2+1)]
            #print(asub)
            dx=getDistanceFromLatLonInKm(lats[0],lons[i1],lats[0],lons[i2])
            if (np.mean(asub)>0):
                area_vals.append(dx)
                anom_vals.append(np.mean(asub))
        else:
            for j in range(0,(ylen-1)):  
                j1=j
                j2=j+1
                #print("i is {:},{:} and j is {:},{:}".format(i1,i2,j1,j2))
                asub=avals[j1:(j2+1),i1:(i2+1)]
                dx1=getDistanceFromLatLonInKm(lats[j1],lons[i1],lats[j1],lons[i2])
                dx2=getDistanceFromLatLonInKm(lats[j2],lons[i1],lats[j2],lons[i2])
                dy=getDistanceFromLatLonInKm(lats[j1],lons[i1],lats[j2],lons[i1])
                if (np.mean(asub)>0):
                    area_vals.append((dx1+dx2)*0.5*dy)
                    anom_vals.append(np.mean(asub))
    #print("resultant array is {:} long".format(len(area_vals)))
    area_wgt_anom = np.array(anom_vals)*np.array(area_vals)
    store_anom['sum_area']=sum(area_vals)
    store_anom['amax'] = amax
    store_anom['AI']=sum(area_wgt_anom)/100000000
    return(store_anom)

def calcBI(field_mask):
    #How many nonzero points are there?
    nnzero=np.count_nonzero(field_mask.values)
    nsample=math.ceil(nnzero/10)
    ind_nonzero = np.nonzero(field_mask.values)
    mask_nonzero = field_mask.values[ind_nonzero]
    ind_sort = np.argsort(mask_nonzero)
    mask_sort = mask_nonzero[ind_sort]
    bot_vals = mask_sort[:nsample]
    top_vals = mask_sort[-nsample:]
    zmin = np.mean(bot_vals)
    zmax = np.mean(top_vals)
    rc = (zmin+zmax)*0.5
    bi = 100.*(zmax/rc-1.)
    return(bi)



parser=argparse.ArgumentParser(description="Provide the list of blobs, list of z500, list of anomaly, and the stats text file")
parser.add_argument("-fm","--filemerged",required=True,action="store")
parser.add_argument("-lb","--listblob",required=True,action="store")
parser.add_argument("-vb","--varblob",default="Z_BLOB")
parser.add_argument("-lz","--listz500",required=True,action="store")
parser.add_argument("-vz","--varz500",default="zg")
parser.add_argument('--is4D',action="store_true")
parser.add_argument("--hpa",action="store_true")
parser.add_argument("--timename",default="time")
parser.add_argument("--levname",default="plev")
parser.add_argument("-ld","--listdevs",required=True,action="store")
parser.add_argument("-vd","--vardevs",default="ADZ")
parser.add_argument("-o","--out",required=True,action="store")

results=parser.parse_args()

fname_out=results.out
fname=results.filemerged
data_blob=results.listblob
data_orig=results.listz500
data_dev=results.listdevs

#xarray of blobfiles
dblob=xa.open_mfdataset(data_blob,use_cftime=True)
dblob=dblob.rename({results.timename:'time',results.latname:'lat',results.lonname:'lon'})
dblob=dblob.sortby(dblob['time'])
var_blob=dblob[results.varblob]
#Get the calendar by not decoding 
temp_open=xa.open_dataset(data_blob[0],decode_times=False)
file_calendar=temp_open[results.timename].attrs['calendar']
#xarray of anom files
ddev=da.open_mfdataset(data_dev,use_cftime=True)
ddev=ddev.rename({results.timename:'time',results.latname:'lat',results.lonname:'lon'})
ddev=ddev.sortby(ddev['time'])
var_devs=ddev[results.vardevs]
#Determine what the input calendar should be
tvar_check=var_devs['time'].values
input_calendar='noleap'
if (isinstance(tvar_check[0],cftime.Datetime360Day)==True):
    input_calendar='360_day'
#xarray of z500 files
#Note that these might have leap days, which need to be removed!!
dorig=xa.open_mfdataset(data_orig,use_cftime=True)
dorig=dorig.rename({results.timename:'time',results.latname:'lat',results.lonname:'lon'})
if (input_calendar!='360_day'):
    dorig=dorig.sel(time=~((dorig['time']dt.month==2) & (dorig['time'].dt.day == 29)))
if (results.is4D==True):
    if (results.levname !='plev'):
        dorig=dorig.rename({results.levname:'plev'})
    pval=50000.
    if (results.hpa==True):
        pval=500.
    var_z500=dorig[results.varz500].sel(plev=pval)
else:
    var_z500=dorig[results.varz500]

#Issue: If it's a standard calendar, it's encoded as Gregorian (we want noleap)
#if it's noleap, it's encoded as noleap
#if it's 360, it's encoded as 360
if (isinstance(tvar_check[0],cftime.DatetimeNoLeap) == False) & (input_calendar=='noleap'):
    #convert the time variable for each of the variables
    t1=var_z500['time'].values
    for t in range(0,len(t1)):
        tstep=t1[t]
        t1[t] = cftime.DatetimeNoLeap(tstep.year,tstep.month,tstep.day,tstep.hour)
    var_z500['time'] = t1
    t2=var_devs['time'].values
    for t in range(0,len(t2)):
        tstep=t2[t]
        t2[t] = cftime.DatetimeNoLeap(tstep.year,tstep.month,tstep.day,tstep.hour)
    var_devs['time'] = t2
    t3=var_blob['time'].values
    for t in range(0,len(t3)):
        tstep=t3[t]
        t3[t] = cftime.DatetimeNoLeap(tstep.year,tstep.month,tstep.day,tstep.hour)
    var_blob['time'] = t3


#Merged file
blob_extents=pd.read_csv(fname)
blob_extents['year']=blob_extents['time'].str.slice(0,4)
blob_extents['month']=blob_extents['time'].str.slice(5,7)
blob_extents['day']=blob_extents['time'].str.slice(8,10)
blob_extents['sec']=blob_extents['time'].str.slice(11,)

tvar_check=var_blob['time'].values
df_indices = pd.DataFrame()
for it,r in blob_extents.iterrows():
    ysnum=int(r['year'])
    msnum=int(r['month'])
    dsnum=int(r['day'])
    ssnum=int(r['sec'])/3600
    if (input_calendar=="360_day"):
        sdate=cftime.Datetime360Day(ysnum,msnum,dsnum,ssnum)
    else:
        sdate=cftime.DatetimeNoLeap(ysnum,msnum,dsnum,ssnum)
        
    xmin=r['minlon']-0.001
    xmax=r['maxlon']+0.001
    ymin=r['minlat']-0.001
    ymax=r['maxlat']+0.001
    if (sdate in tvar_check):
        #Need to deal with the case where there's the periodic boundary
        if (xmin>xmax):
            b1=var_blob.sel(time=sdate,lat=slice(ymin,ymax),lon=slice(xmin,360))
            b2=var_blob.sel(time=sdate,lat=slice(ymin,ymax),lon=slice(0,xmax))
            bslice=xa.concat([b1,b2],dim='lon')
            d1=var_devs.sel(time=sdate,lat=slice(ymin,ymax),lon=slice(xmin,360))
            d2=var_devs.sel(time=sdate,lat=slice(ymin,ymax),lon=slice(0,xmax))
            dslice=xa.concat([d1,d2],dim='lon')
            o1=var_z500.sel(time=sdate,lat=slice(ymin,ymax),lon=slice(xmin,360))
            o2=var_z500.sel(time=sdate,lat=slice(ymin,ymax),lon=slice(0,xmax))
            oslice = xa.concat([o1,o2],dim='lon')
        else:
            oslice = var_z500.sel(time=sdate,lat=slice(ymin,ymax),lon=slice(xmin,xmax))
            bslice=var_blob.sel(time=sdate,lat=slice(ymin,ymax),lon=slice(xmin,xmax))
            dslice=var_devs.sel(time=sdate,lat=slice(ymin,ymax),lon=slice(xmin,xmax))
        bmask=bslice*dslice
        omask=bslice*oslice
        rnew = calcAI(bmask,r)
        rnew['BI']=calcBI(omask)
        df_indices=df_indices.append(rnew,ignore_index=True)


df_indices.to_csv(fname_out,index=False,na_rep="_")  

