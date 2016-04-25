from pandas import DataFrame
import numpy

#Pandas indexing:
#  df[col] = selects by column name
#  df.loc[ser] = selects by row name
#  df.iloc[ser] = selects by row index
#  


class BlobObject:
  def __init__(self, id, colnames, it, nlines,lst):
    self.id = id
    self.colnames = colnames
    self.data = DataFrame(index=range(it,it+nlines), columns=colnames)
    for l in lst:
      line = map(float,l.strip().split())
      #first index is the row number index
      rn=int(line[0])
      info=line[1:]
      self.data.set_value(rn,colnames,info)


    

f = open("/global/cscratch1/sd/marielp/2xCO2/data/y_16_stats","r")
lines=f.readlines()
print("Read file.")
#Search for all instances of "Blob"

BlobIndices=[]

for i,j in enumerate(lines):
  if "Blob" in j:
    BlobIndices+=[i]
print("Saved indices")

colnames = ['minlat','maxlat','minlon','maxlon','centlat','centlon','area']

print("entering loop")
#Now, for each blob,create a BlobObject
for n in range(0,3):
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
  Blob.data.boxplot(column='minlat')


