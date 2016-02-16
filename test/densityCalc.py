from netCDF4 import Dataset
import numpy as np
import sys

#nc_filename = sys.argv[1]
nc_filename = '/media/mariellep/my_hd/ERA_data/ERA_2013_MAM_5day_blobs.nc'
nc_out = nc_filename.replace("blobs.nc","density.nc")
nc=Dataset(nc_filename,mode='r')

blobs=nc.variables['INT_ADIPVtag'][:]
lon=nc.variables['lon'][:]
lat=nc.variables['lat'][:]
t=nc.variables['time'][:]
lat_len = len(lat)
lon_len = len(lon)
#nc.close()

div = np.divide(blobs,blobs)

dens=np.sum(blobs,axis=0)
denom = float(len(t))
dens=dens/denom

dens_out = Dataset(nc_out,mode='w')
dens_out.createDimension('lat',lat_len)
dens_lat = dens_out.createVariable('lat',nc.variables['lat'].dtype,('lat',))
for a in nc.variables['lat'].ncattrs():
  dens_lat.setncattr(a,nc.variables['lat'].getncattr(a))
dens_out.variables['lat'][:] = nc.variables['lat'][:]

dens_out.createDimension('lon',lon_len)
dens_lon=dens_out.createVariable('lon',nc.variables['lon'].dtype,('lon',))
for a in nc.variables['lon'].ncattrs():
  dens_lon.setncattr(a,nc.variables['lon'].getncattr(a))
dens_out.variables['lon'][:] = nc.variables['lon'][:]

dens_write = dens_out.createVariable('dens','f8',('lat','lon'))
dens_write[:] = dens[:]

nc.close()
dens_out.close()
