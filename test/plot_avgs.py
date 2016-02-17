from netCDF4 import Dataset
import numpy as np
import sys
import matplotlib.pyplot as plt

def month_len(month_num):
  if month_num==4 or month_num==6 or month_num==9 or month_num==11:
    num_days = 30
  elif month_num ==2:
    num_days = 28
  else:
    num_days = 31
  return num_days


latindex = 30
lonindex = 0

davgname = '/media/mariellep/ERA_data/ERA_avg/ERA_1981_2005_dailyavg.nc'
davg = Dataset(davgname,'r')


