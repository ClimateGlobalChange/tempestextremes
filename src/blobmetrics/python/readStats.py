import argparse
import pandas as pd
pd.set_option('mode.chained_assignment',None)

#THIS READS IN BLOBSTATS DATA AND OUTPUTS IT IN CSV FORMAT
parser=argparse.ArgumentParser(description="Provide a text file list of BlobStats outputs")
parser.add_argument("-f","--filelist",required=True,action="store")
parser.add_argument("-o","--out",required=True,action="store")
parser.add_argument("-s","--season",required=True,action="store")
parser.add_argument("-r","--region",required=True,action="store")
parser.add_argument("-d","--dataset",required=True,action="store")
#parser.add_argument("-c","--calendar",required=True,action="store")
results=parser.parse_args()
flist=open(results.filelist).read().splitlines()
s=results.season
r=results.region
d=results.dataset
fname_out=results.out

#Get the headers from the first file
colnames=open(flist[0]).readlines()[0].strip().split(",")
colnames_full=["var","bnum","season","region"]
colnames_full.extend(colnames)
colnames_full.append("fname")

#Output dataframe
dat_total=pd.DataFrame(columns=colnames_full)
for str_file in flist:
    dat=open(str_file).readlines()
    nl=len(dat)
    if(nl>1):
        for l in range(1,nl):
            lin=dat[l].strip().split("\t")
            if "Blob" in lin[0]:
                bnum=lin[0].split()[1]
            else:
                dict_line={"var":d,"bnum":bnum,"time":lin[0],"season":s,"region":r,"fname":str_file}
#                dict_line={"var":d,"bnum":bnum,"time":lin[0],"season":s,"region":r,"calendar":results.calendar,"fname":str_file}
                nc=1
                for v in colnames[1:len(colnames)]:
                    dict_line[v]=float(lin[nc])
                    nc+=1
                dat_total=dat_total.append(dict_line,ignore_index=True)

dat_total['area_km'] = dat_total['area']*(6371.*6371.)
colnames_dat=dat_total.columns.tolist()
colnames_dat[len(colnames_dat)-2]='area_km'
colnames_dat[len(colnames_dat)-1]='fname'
#reorder the columns (it was bugging me)
dat_total=dat_total[colnames_dat]
dat_total.to_csv(fname_out,index=False,na_rep='_')
print("wrote {:}".format(fname_out))
