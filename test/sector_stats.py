#!/bin/python

import numpy as np
from netCDF4 import Dataset
import sys
import matplotlib.pyplot as plt

def get_index(arr,n1,tol):
  y1=np.where(abs(n1-arr)<tol)
  n1=y1[0][0]
  return(n1)

def get_arr(var,sec,hemi):
  if hemi=="NH":
    l1=get_index(lat,15.,0.5)
    l2=get_index(lat,75.,0.5)
  
    if sec=="ATL":
      n1=get_index(lon,280.,.5)
      n2=get_index(lon,40.,.5)    
    elif sec=="PAC":
      n1=get_index(lon,140.,0.5)
      n2=get_index(lon,260.,0.5)
    elif sec=="CONT1":
      n1=get_index(lon,261.,0.5) 
      n2=get_index(lon,279.,0.5)
    elif sec=="CONT2":
      n1=get_index(lon,41.,0.5)
      n2=get_index(lon,139.,0.5)

  elif hemi=="SH":
    l1=get_index(lat,-75.,0.5)
    l2=get_index(lat,-15.,0.5)
   
    if sec=="ATL":
      n1=get_index(lon,300.,0.5)
      n2=get_index(lon,30.,0.5)
    elif sec=="PAC":
      n1=get_index(lon,130.,0.5)
      n2=get_index(lon,300.,0.5)
    elif sec=="IO":
      n1=get_index(lon,31.,0.5)
      n2=get_index(lon,129.,0.5)

  if sec=="ATL":
    ATL1=var[l1:l2,n1:]
    ATL2=var[l1:l2,:n2]
    outarr=np.concatenate((ATL1,ATL2),axis=1)
  else:
   outarr=var[l1:l2,n1:n2]

  return(outarr)


############################################################
#nc1name=sys.argv[1]
#nc2name=sys.argv[2]

nc1name="/global/cscratch1/sd/marielp/climo/blobs/JJA_avg_dens_climo.nc"
nc2name="/global/cscratch1/sd/marielp/2xCO2/blobs/JJA_avg_dens_2xCO2.nc"

nc1=Dataset(nc1name,"r")
nc2=Dataset(nc2name,"r")

nc1dens=nc1.variables['dens'][:]
nc2dens=nc2.variables['dens'][:]

lon=nc1.variables['lon'][:]
lat=nc1.variables['lat'][:]

#ATL: -80. (280.) to 40. (NH), -60. (300.) to 30. (SH)
#PAC: 140. to -100. (260.) (NH), 130. to -60. (300.)(SH)
#Continental: -100. to -80. and 40. to 140. (NH)
# Indian Ocean: 30. to 130.

NH_ATL1=get_arr(nc1dens,"ATL","NH")
NH_ATL2=get_arr(nc2dens,"ATL","NH")


