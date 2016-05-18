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

#Note: current code starts with 0016-01-01 and ends with 0016-12-31
# For seasonal or other years, might have some issues with first December file
# (starts on 12-07) or January file (0018 starts with 01-06)
        
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

#plt.figure(2)
#plt.figure(3)
print("entering loop")
NHVec = []
SHVec = []

minarea = 10e20
maxarea = 0.
maxt=0.
#Now, for each blob,create a BlobObject
for n in range(0,len(BlobIndices)):
#for n in range(1,5):
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
    
  if (Blob.data.t_since_init.max()>maxt):
    maxt=Blob.data.t_since_init.max()
    
  if Blob.hemi == 'NH':
    NHVec += [Blob]
  else:
    SHVec += [Blob]
  #plt.figure(1)
  #Blob.map_latlon('centlon','centlat',mp)
  #plt.figure(2)
  #Blob.plot_mvmt()
  #plt.figure(3)
  #Blob.plot_area()


#plt.figure(figsize=(24,12))
plt.figure(1)
#setup Basemap
mp=Basemap(projection='cyl',llcrnrlat=-90.,urcrnrlat=90.,
            llcrnrlon=0,urcrnrlon=360.,resolution='c')
mp.drawcoastlines()
mp.drawparallels(numpy.arange(-90.,91.,30.),labels=[0,0,0,1])
mp.drawmeridians(numpy.arange(0.,360.,60.),labels=[0,0,0,1])

#plt.figure(2)
#n=44 (Blob 85) and n=51 (Blob 99) have huge jumps in longitude!
n=0
for b in NHVec:
  print "n is %d"%n
  print b.name
  plt.figure(1)
  map_t_latlon(b.data,maxt,mp)
  n+=1
  #plt.plot(b.data.index.values,b.data.centlon.values)


#Some histograms
#plt.figure(4)
#plt.figure(3)
#plt.subplot(211)
#max_area = [b.data.area.values.max() for b in NHVec]
#dist = [b.dist for b in NHVec]
#plt.scatter(dist,max_area)
#plt.hist(init_lon)
#plt.subplot(212)
#max_area = [b.data.area.values.max() for b in SHVec]
#dist = [b.dist for b in SHVec]
#plt.scatter(dist,max_area)

#plt.figure(5)
#plt.subplot(211)
#ndays = [b.ndays for b in NHVec]
#plt.hist(ndays)
#plt.subplot(212)
#ndays = [b.ndays for b in SHVec]
#plt.hist(ndays)
