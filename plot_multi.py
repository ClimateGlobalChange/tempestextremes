from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import sys


#import relevant variables
nc_start = sys.argv[1]
nc_year1=int(sys.argv[2])
nc_year2=int(sys.argv[3])

nc=Dataset(nc_start,mode='r')

lon=nc.variables['lon'][:]
lat=nc.variables['lat'][:]
dens=nc.variables['dens'][:]

nc.close()

#range of years in average
nc_year1+=1
nc_year2+=1

a=1.
nc_namesplit=nc_start.split("_")
for n in range(nc_year1,nc_year2):
  nc_copy = nc_namesplit
  for s in range(0,len(nc_copy)):
    if nc_copy[s].isdigit():
      nc_copy[s] = str(n)
  nc_name = "_".join(nc_copy)
  print nc_name
  nc2=Dataset(nc_name,mode='r')
  d2=nc2.variables['dens'][:]
  nc2.close()
  dens +=d2
  a+=1.

dens=dens/a

print nc_namesplit
for s in range (0,len(nc_namesplit)):
  if nc_namesplit[s].isdigit():
    nc_namesplit[s] = "%d_%d"%(nc_year1-1,nc_year2-1)
nc_split2 = "_".join(nc_namesplit)
nc_split3 = nc_split2.split(".")

nc_out = nc_split3[0] + "_plot.png"
lat_ticks=np.arange(-90,90,30)
lon_ticks=np.arange(0,360,30)

#create map
plt.figure(figsize=(24,12))
m = Basemap(llcrnrlat=-90.,urcrnrlat=90.,
            llcrnrlon=0.,urcrnrlon=359.,
            resolution='c')
         
m.drawcoastlines()
m.drawmapboundary()
m.drawparallels(lat_ticks,labels=[0,0,0,1])
m.drawmeridians(lon_ticks,labels=[0,0,0,1])

x,y = m(*np.meshgrid(lon,lat))
#draw filled contours
pv=m.contourf(x,y,dens,levels=np.arange(0,0.25,.025),extend='both')
#add colorbar
cb=plt.colorbar(pv)
plt.title("PV* density for %s %d to %d"%(sys.argv[4],nc_year1-1,nc_year2-1))
plt.savefig(nc_out,bbox_inches='tight')
