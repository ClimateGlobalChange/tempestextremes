#!/bin/python
import cdms2
import sys
import numpy
import MV2

# FLAGS
cdms2.setNetcdfShuffleFlag(0)
cdms2.setNetcdfDeflateFlag(0)
cdms2.setNetcdfDeflateLevelFlag(0)
file_in = sys.argv[1]
#varname = sys.argv[2]
file_open=cdms2.open(file_in,'r')
#file_open=cdms2.open("ERA_2013_blobs.nc",'r')
varname = "PV_ADEV_INTtag"
tag=file_open.variables[varname]
time_axis = file_open.axes['time']
tag_arr=numpy.zeros(tag.shape)
tag_arr[:]=tag[:]
with numpy.errstate(divide='ignore'):
  num=tag_arr
  den=tag_arr
  result=num/den
  result[den==0] = 0


int_array = MV2.array(result)
int_array.id='BIN'
int_array.setAxis(0,time_axis)
int_array.setAxis(1,file_open.axes['lat'])
int_array.setAxis(2,file_open.axes['lon'])


num_steps = time_axis.shape[0]
tag_avg = numpy.sum(result,0)/num_steps
avg_array = MV2.array(tag_avg)
avg_array.id='DENS'
avg_array.setAxis(0,file_open.axes['lat'])
avg_array.setAxis(1,file_open.axes['lon'])


fname_out = file_in.replace(".nc","_dens.nc")
#fname_out='debug.nc'
file_out = cdms2.open(fname_out, "w")
for attributes in file_open.listglobal():
	setattr(file_out, attributes, getattr(file_open, attributes))
file_out.write(avg_array)
file_out.write(int_array)
file_out.close()

#NCRA command: ncra -d time,0,,1 [file_out] [new file out]
