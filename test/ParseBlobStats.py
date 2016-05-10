from pandas import DataFrame
import pandas
import numpy
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.collections import LineCollection
from mpl_toolkits.basemap import Basemap, addcyclic
import math
from ParseBlobFuncs import *


pandas.options.mode.chained_assignment = None
#Pandas indexing:
#  df[col] = selects by column name
#  df.loc[row] = selects by row name
#  df.iloc[ind] = selects by row index#  


#Bounds of colormap and the colors themselves
#cnorm = mpl.colors.Normalize(vmin, vmax)
#create a cmin and cmax default
# mpl.cm.ScalarMappable(cmap=pylab.get_cmap('rainbow'))
# c=[scalapMap.to_rgba(variable)]


#Note: current code starts with 0016-01-01 and ends with 0016-12-31
# For seasonal or other years, might have some issues with first December file
# (starts on 12-07) or January file (0018 starts with 01-06)
def BlobStats(BlobVec, plot_type, start_lon, end_lon):
  start_lon = 0.
  end_lon = 360.
  
  per_wrap = False
  
  if start_lon<0.:
    start_lon += 360.
  if end_lon<0.:
    end_lon += 360.
  if start_lon > end_lon:
    per_wrap = True
    
def plot_all(nvec=None,svec=None,bmap=None):
    if ((nvec != None) and (svec != None)):
      vec = nvec + svec
    elif (svec == None):
      vec = nvec
    else:
      vec = svec
    

      
    
###############################################
#f = open("/global/cscratch1/sd/marielp/2xCO2/data/y_16_stats","r")
f = open("y_03_stats_2xCO2","r")
lines=f.readlines()
print("Read file.")
#Search for all instances of "Blob"
BlobIndices=[]

for i,j in enumerate(lines):
  if "Blob" in j:
    BlobIndices+=[i]
print("Saved indices")

colnames = ['minlat','maxlat','minlon','maxlon','centlat','centlon','area']

plt.figure(1)
#setup Basemap
mp=Basemap(projection='cyl',llcrnrlat=-90.,urcrnrlat=90.,
            llcrnrlon=0,urcrnrlon=360.,resolution='c')
mp.drawcoastlines()
mp.drawparallels(numpy.arange(-90.,91.,30.))
mp.drawmeridians(numpy.arange(0.,360.,60.))

#plt.figure(2)
#plt.figure(3)
print("entering loop")
NHVec = []
SHVec = []

minarea = 10e20
maxarea = 0.
#Now, for each blob,create a BlobObject
for n in range(0,len(BlobIndices)):
#for n in range(12,30):
  if (n < len(BlobIndices)-1):
    subarr = lines[BlobIndices[n]:BlobIndices[n+1]]
  else:
    subarr = lines[BlobIndices[n]:]

  BlobTitle=subarr[0].strip().split()
  firstline=map(float,subarr[1].strip().split())

  i = int(firstline[0])
  title = BlobTitle[1]
  nl = int(filter(str.isdigit,BlobTitle[2]))

  Blob=BlobObject(title,colnames,i,nl,subarr[1:])
  if (Blob.data.area.min()<minarea):
    minarea=Blob.data.area.min()
  if (Blob.data.area.max()>maxarea):
    maxarea=Blob.data.area.max()
#  if Blob.hemi == 'NH':
#    NHVec += [Blob]
#  else:
#    SHVec += [Blob]
  plt.figure(1)
  Blob.map_latlon('centlon','centlat',mp)
  #plt.figure(2)
  #Blob.plot_latlon('centlon','centlat')
  #plt.figure(3)
  #Blob.plot_area()

#Some histograms
#plt.figure(4)
#plt.subplot(211)
#init_lon = [b.data.centlon.iloc[0] for b in NHVec]
#plt.hist(init_lon)
#plt.subplot(212)
#init_lon = [b.data.centlon.iloc[0] for b in SHVec]
#plt.hist(init_lon)

#plt.figure(5)
#plt.subplot(211)
#ndays = [b.ndays for b in NHVec]
#plt.hist(ndays)
#plt.subplot(212)
#ndays = [b.ndays for b in SHVec]
#plt.hist(ndays)
