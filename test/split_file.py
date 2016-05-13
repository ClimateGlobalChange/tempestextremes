#!/bin/python

from netCDF4 import Dataset
import sys

fname= sys.arvg[1]
fout1 = sys.argv[2]
fout2 = sys.argv[3]


f=Dataset(fname,"r")
tdim = f.dimensions['time']
tvar = f.variables['time']

tsplitlen = len(tvar)/2
if (len(tvar)%2 != 0):
  print "Axis has odd number of time steps"



