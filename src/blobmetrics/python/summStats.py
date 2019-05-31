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


parser=argparse.ArgumentParser(description="Parse the merged/index file and summarize")
parser.add_argument("-f","--filein",required=True,action="store")
parser.add_argument("-o","--out",required=True,action="store")
parser.add_argument("-c","--calendar",required=True,action="store")
results=parser.parse_args()

calendar=results.calendar
fname_out=results.out
dat_in=read.csv(results.filein,na_filter=False)
dat_summ=pd.DataFrame()

for f in sorted(list(set(dat_in['fname']))):
    dat_sub = dat_in[dat_in['fname'].str.match(f)]
    bnum_unique=set(dat_sub['bnum'])
    for b in sorted(list(bnum_unique)):
        dat_bsub=dat_sub.loc[dat_sub['bnum']==b]
        dat_bsub = dat_bsub.sort_values('time')
        sline=dat_bsub.iloc[0]
        eline=dat_bsub.iloc[len(dat_bsub)-1]
        stime=sline['time']
        etime=eline['time']
        ysnum=int(stime[:4])
        msnum=int(stime[5:7])
        dsnum=int(stime[8:10])
        ssnum=int(stime[11:])/3600
        yenum=int(etime[:4])
        menum=int(etime[5:7])
        denum=int(etime[8:10])
        senum=int(etime[11:])/3600
        if calendar=="360_day":
            sdate=cftime.Datetime360Day(ysnum,msnum,dsnum,ssnum)
            edate=cftime.Datetime360Day(yenum,menum,denum,senum)
        else:
            sdate=cftime.DatetimeNoLeap(ysnum,msnum,dsnum,ssnum)
            edate=cftime.DatetimeNoLeap(yenum,menum,denum,senum)
        td=edate-sdate
        num_days = td.days +1
        num_hrs = td.seconds/(3600*24)
        num_tot_days = num_days + num_hrs 
        avg_clat=(sline['centlat']+eline['centlat'])/2.
        dict_bsub={"var":sline['var'],"bnum":b,"region":sline['region'],"startdate":sline['time'],"enddate":eline['time'],
                    "duration_days":num_tot_days,"mean_centlat":dat_bsub['centlat'].mean(),
                    "start_centlat":sline['centlat'],"start_centlon":sline['centlon'],
                    "end_centlat":eline['centlat'],"end_centlon":eline['centlon'],
                    "mean_centlat":dat_bsub['centlat'].mean(),
                    "dist_km":getDistanceFromLatLonInKm(sline['centlat'],sline['centlon'],
                                                        eline['centlat'],eline['centlon']),
                    "zonal_dist_km":getDistanceFromLatLonInKm(avg_clat,sline['centlon'],
                                                            avg_clat,eline['centlon']),
                    "min_area_km":dat_bsub['area_km'].min(),"max_area_km":dat_bsub['area_km'].max(),
                    "mean_area_km":dat_bsub['area_km'].mean(),"min_AI":dat_bsub['AI'].min(),
                    "max_AI":dat_bsub['AI'].max(),"mean_AI":dat_bsub['AI'].mean(),
                    "min_BI":dat_bsub['BI'].min(),"max_BI":dat_bsub['BI'].max(),"mean_BI":dat_bsub['BI'].mean()}
        dict_bsub['zonal_speed_km']=dict_bsub['dist_km']/(dict_bsub['duration_days']*24.)
        dict_bsub['fname']=sline['fname']
        if (len(dat_bsub[dat_bsub['merged'].str.contains('YES')])>0):
            dict_bsub['merged'] = 'YES'
        else:
            dict_bsub['merged'] = 'NO'
        dat_summ=dat_summ.append(dict_bsub,ignore_index=True)
#Write to file
colnames_summ=dict_bsub.keys()
dat_summ=dat_summ[colnames_summ]

dat_summ.to_csv(fname_out,index=False,na_rep="_")


