#This script re-formats the ERA-Interim data so it can be properly read by VisIt
import sys
import MV2
import cdms2
from array import array

cdms2.setNetcdfShuffleFlag(0)
cdms2.setNetcdfDeflateFlag(0)
cdms2.setNetcdfDeflateLevelFlag(0)

file_in = cdms2.open(sys.argv[1], 'r')

time_axis = file_in.axes['initial_time0_hours']
time_axis.id = 'time'
time_units = time_axis.units
DefaultCalendar=time_axis.getCalendar()
time_axis.designateTime()

plev_axis = file_in.axes['lv_ISBL1']
plev_axis.id = 'lev'
plev_axis.designateLevel()

lat_axis = file_in.axes['g0_lat_2']
lat_axis.id = 'lat'
lat_axis.designateLatitude()

lon_axis = file_in.axes['g0_lon_3']
lon_axis.id = 'lon'
lon_axis.designateLongitude()

outname = sys.argv[1].replace(".nc", "_mod.nc")
file_out = cdms2.open(outname, 'w')

for attributes in file_in.listglobal():
  setattr(file_out, attributes, getattr(file_in, attributes))

t_var = file_in.variables['T_GDS0_ISBL']
t_var.id = 'T'

u_var = file_in.variables['U_GDS0_ISBL']
u_var.id = 'U'

v_var = file_in.variables['V_GDS0_ISBL']
v_var.id = 'V'

#vort_var = file_in.variables['VO_GDS0_ISBL']
#vort_var.id = 'VO'

#pv_var = file_in.variables['PV_GDS0_ISBL']
#pv_var.id = 'PV'

#z_var = file_in.variables['Z_GDS0_ISBL']
#z_var.id = 'Z'

#var_dict = {'T':t_var, 'U':u_var, 'V':v_var, 'VO':vort_var, 'PV':pv_var, 'Z':z_var}
var_dict = {'T':t_var,'U':u_var,'V':v_var}
var_keys = var_dict.keys()
var_len = len(var_keys)

var_shape = (time_axis.shape[0], plev_axis.shape[0], lat_axis.shape[0], lon_axis.shape[0])
for i in range(0, var_len):
  var_name = var_keys[i]
  var_call = MV2.array(var_dict[var_name])
  if var_call.shape == var_shape:
    var_call.setAxis(0, time_axis)
    var_call.setAxis(1, plev_axis)
    var_call.setAxis(2, lat_axis)
    var_call.setAxis(3, lon_axis)
    file_out.write(var_call)

file_out.close()
file_in.close()
sys.exit(0)
