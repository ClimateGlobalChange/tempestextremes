from pandas import DataFrame
import pandas
import numpy
import pylab
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from mpl_toolkits.basemap import Basemap, addcyclic
import math

def great_circle_dist(lat1,lon1,lat2,lon2):
  r = 6.371e6
  
  lat1=math.radians(lat1)
  lon1=math.radians(lon1)
  lat2=math.radians(lat2)
  lon2=math.radians(lon2)
  dlat2=0.5*(lat1-lat2)
  dlon2=0.5*(lon1-lon2)
  
  a = math.sin(dlat2)*math.sin(dlat2) + \
            math.cos(lat1)*math.cos(lat2)*\
            math.sin(dlon2)*math.sin(dlon2)
  c = 2. * math.atan2(math.sqrt(a),math.sqrt(1.-a))
  d = c * r
  return(d)
  
#ATL: -80. (280.) to 40. (NH), -60. (300.) to 30. (SH)
#PAC: 140. to -100. (260.) (NH), 130. to -60. (300.)(SH)
#Continental: -100. to -80. and 40. to 140. (NH)
# Indian Ocean: 30. to 130.

#Sectors are defined on 0.-360. range since this is what the axis is defined on!
def sector(hemi,lon0):
  if hemi == "NH":
    if ((lon0 < 40.) or (lon0 > 280.)):
      sec = "ATL"
    elif ((lon0 > 140.) and (lon0 < 260.)):
      sec = "PAC"
    else:
      sec = "CONT"
      
  elif hemi == "SH":
    if ((lon0 < 30.) or (lon0 > 300.)):
      sec = "ATL"
    elif ((lon0 > 130.) and (lon0 < 300.)):
      sec = "PAC"
    else:
      sec = "IO"
    #print("for lon0 %f, sector is %s"%(lon0,sec))
  return(sec)
    


#Function that will generate color values for each point (on a scale of 0-1)
def color_split_num(var, minvar, maxvar, n=10):
  ncols = []

  bins=numpy.linspace(minvar,maxvar,n)
  col_idx = numpy.linspace(0,1,n)
  vals = var.astype('float64')

  #This will return the index of the bin value for each item
  d = numpy.digitize(vals, bins)
  for i in range (0,len(var)):
    ncols += [col_idx[d[i]-1]]

  return(ncols)
  
#Function that will generate color values for each unique string (on a scale of 0-1)
#uniq_var=set(var)
def color_split_str(var, uniq_var):
  ncols = []

  n_uniq = len(uniq_var)
  list_var = list(uniq_var)
  col_idx = numpy.linspace(0,1,n_uniq)

  for i in range(0,len(list_var)):
    for j in range(0,len(var)):
      if list_var[i] == var[j]:
        ncols+= [col_idx[i]]
    
  return(ncols)

def plot_mcline(var1,var2,cvar,vartype,xmin,xmax,ymin,ymax,cmin=0,cmax=1,n=10,c='viridis'):
                
  points = numpy.array([var1.astype('float64'),\
             var2.astype('float64')]).T.reshape(-1, 1, 2)
  segments = numpy.concatenate([points[:-1], points[1:]], axis=1) 
  if (vartype == 'STR'):
    uniq_var=set(cvar)
    cols = color_split_str(cvar,uniq_var)
  elif (vartype == 'NUM'):
    cols = color_split_num(cvar,cmin,cmax,n)
   
  dummy_cols = numpy.asarray(cols)
  
  lc = LineCollection(segments, cmap=cm.get_cmap(c))
  lc.set_array(dummy_cols)
  lc.set_linewidth(1.5)
  plt.gca().add_collection(lc)
  plt.xlim(xmin,xmax)
  plt.ylim(ymin, ymax)

def map_t_latlon(d,tmax,m):
  #vals=color_split_num(d['t_since_init'].values,cmin,cmax)
  dat=d[['centlon','centlat','t_since_init']]
  #dat['col']=vals
  #print(dat)
  
  #Pick out areas where there's a periodic boundary condition
  #note that these are not index values!!! 
  per_values = dat[abs(dat.centlon.diff()) > 100. ].index.tolist()
  print per_values
  #x,y=m(dat.centlon.values.astype('float64'),dat.centlat.values.astype('float64'))
  if (len(per_values) == 0):
    print "No periodic condition"
    print dat.centlon.values
    x,y=m(dat.centlon.values.astype('float64'),dat.centlat.values.astype('float64'))

    plot_mcline(x,y,dat['t_since_init'].values,'NUM',\
    0.,360.,-90.,90.,0,int(tmax),n=int(tmax))

  else:
    startloc=dat.index[0]
    for n in range(0,len(per_values)):  
      endloc=per_values[n]-1
      dat_sub = dat.loc[startloc:endloc]
      print dat_sub.centlon.values
      startloc=per_values[n]

    dat_sub=dat.loc[startloc:]
    print dat_sub.centlon.values
      #x,y=m(dat_sub.centlon.values, dat_sub.centlat.values)
      #plot_mcline(x,y,dat_sub['t_since_init'].values,'NUM',\
       #0.,360.,-90.,90.,0,int(tmax),n=int(tmax))
      

#Plotting functions
def plot_latlon(dat,latname,lonname):
  clat=dat[latname]
  clon=dat[lonname]
  #Pick out areas where there's a periodic boundary condition
  #note that these are not index values!!! 
  per_values = dat[abs(clon.diff()) > 300 ].index.tolist()
  if (len(per_values) == 0):
    plt.plot(clon.values,clat.values)
  else:
    for n in range(0,len(per_values)):
      if (n < len(per_values)-1):
        clat_sub = clat.loc[per_values[n]:per_values[n+1]-1]
        clon_sub = clon.loc[per_values[n]:per_values[n+1]-1]
      else:
        clat_sub = clat.loc[per_values[n]:]
        clon_sub = clon.loc[per_values[n]:]
      
      plt.plot(clon_sub.values,clat_sub.values)
      plt.ylabel(latname)
      plt.xlabel(lonname)  
  
  
class BlobObject:
  #initialize object:
  #name of blob
  #column names
  #starting index number
  #number of rows
  #array with data
  def __init__(self, name, colnames, it, nlines,lst):
    self.name = name
    self.colnames = colnames
    self.data = DataFrame(index=range(it,it+nlines), columns=colnames)
    for l in lst:
      line = map(float,l.strip().split())
      #first item in list corresponds to row number index (i.e. the time step)
      rn=int(line[0])
      #The rest of the list is everything else
      info=line[1:]
      self.data.set_value(rn,colnames,info)

    self.averages = self.data.mean()
    self.std_dev = self.data.std()
    
    #Calculate number of days from init
    itime=self.data.index.values[0]
    self.data['t_since_init']=(self.data.index.values-itime)*0.125
    #Calculate distance of center from start to finish    
    clat0=self.data['centlat'].iloc[0]
    clon0=self.data['centlon'].iloc[0]
    clat1=self.data['centlat'].iloc[nlines-1]
    clon1=self.data['centlon'].iloc[nlines-1]
    self.dist=great_circle_dist(clat0,clon0,clat1,clon1)
    if clat0 < 0:
      self.hemi = "SH"
    else:
      self.hemi = "NH"
      
    self.data['sector'] = [sector(self.hemi, l) for l in self.data.centlon]
    self.ndays = float(len(self.data.index))/8.
    
  def plot_mvmt(self):
    aarr = self.data[['centlon','sector','area']]  
    points = numpy.array([aarr.index.values.astype('float64'),\
                aarr.centlon.values.astype('float64')]).T.reshape(-1, 1, 2)
    segments = numpy.concatenate([points[:-1], points[1:]], axis=1) 
    cols = color_split(aarr.area.values,'NUM')
    dummy_cols = numpy.asarray(cols)
    cmin=aarr.centlon.values.min()
    cmax=aarr.centlon.values.max()
    xmin=aarr.index.values.min()
    xmax=aarr.index.values.max()
    print(cmin, cmax)
    print(xmin,xmax)
    #cnorm=mpl.colors.Normalize(cmin,cmax)
    lc = LineCollection(segments, cmap=cm.get_cmap('Dark2_r'))#,\
    #        norm=cnorm )
    lc.set_array(dummy_cols)
    lc.set_linewidth(3)
    #lc.set_array(aarr.values[:-1])
    plt.gca().add_collection(lc)
    plt.xlim(xmin,xmax)
    plt.ylim(cmin-0.01, cmax+0.01)

    plt.ylabel('Longitude')
    plt.xlabel('Time index')


      

