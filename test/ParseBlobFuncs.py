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
def color_split(var, vartype, n=10):
  ncols = []
  if vartype == 'STR':
    uniq_var = set(var)
    n_uniq = len(uniq_var)
    list_var = list(uniq_var)
    col_idx = numpy.linspace(0,1,n_uniq)

    for i in range(0,len(list_var)):
      for j in range(0,len(var)):
        if list_var[i] == var[j]:
          ncols+= [col_idx[i]]
  elif vartype == 'NUM':
    minvar = min(var)
    maxvar = max(var)
    bins=numpy.linspace(minvar,maxvar,n)
    #print(bins)
    col_idx = numpy.linspace(0,1,n)
    #print(col_idx)
    vals = var.astype('float64')
    #print(vals)
    #This will return the index of the bin value for each item
    d = numpy.digitize(vals, bins)
    #print(d)
    for i in range (0,len(var)):
      ncols += [col_idx[d[i]-1]]
    
  #str_ncols = [str(x) for x in ncols]
  #str_ncols = str_ncols[:-1]
  #return(str_ncols)
#  return(col_idx)
  return(ncols)
  
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
    
  #Plotting functions
  def plot_latlon(self,lonname,latname):
    clat=self.data[latname]
    clon=self.data[lonname]
    #Pick out areas where there's a periodic boundary condition
    #note that these are not index values!!! 
    per_values = self.data[abs(clon.diff()) > 300 ].index.tolist()
    plt.axis
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

  def map_latlon(self,lonname,latname,m):
    dat=self.data[[lonname,latname,'area']]
    vals = color_split(dat.area.values, 'NUM')
    
    #Pick out areas where there's a periodic boundary condition
    #note that these are not index values!!! 
    per_values = dat[abs(dat[lonname].diff()) > 300 ].index.tolist()
    if (len(per_values) == 0):
      plt.plot(dat[lonname].values,dat[latname].values)
    else:
      for n in range(0,len(per_values)):
        if (n < len(per_values)-1):
          dat_sub = dat.loc[per_values[n]:per_values[n+1]-1]
        else:
          dat_sub = dat.loc[per_values[n]:]
        x,y=m(dat_sub[lonname].values, dat_sub[latname].values)
        plt.plot(x,y)
        plt.ylabel(latname)
        plt.xlabel(lonname)
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
    plt.colorbar(dummy_cols)


def plot_func(vartype,var1,var2=None,cvar=None,cm='viridis',bmap=None):    
    plt.figure(figsize=(24,12))
    if (vartype=='HIST'):
      plt.hist(var1)
    elif(vartype=='SCTR'):
      if (cvar == None):
        plt.scatter(var1,var2)
      else:
        plt.scatter(var1,var2,c=cvar,cmap=cm)
        plt.colorbar()
      

