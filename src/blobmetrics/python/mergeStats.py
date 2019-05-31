import argparse
import pandas as pd
pd.set_option('mode.chained_assignment',None)

#Look at when the data is merged-- are there any blobs that should be split?
def lon_convert_360_180(lon):
    distFrom180=lon-180
    retlon=-999
    if (distFrom180<0):
        retlon=lon
    else:
        retlon=distFrom180-180.
    return(retlon)


parser=argparse.ArgumentParser(description="Provide the BlobStats outputs from both StitchBlobs and DetectBlobs")
parser.add_argument("-fs","--filestitch",required=True,action="store")
parser.add_argument("-fd","--filedetect",required=True,action="store")
parser.add_argument("-o","--out",required=True,action="store")

results=parser.parse_args()

f_stitch=results.filestitch
f_detect=results.filedetect
fname_out=results.out

df_stitch=pd.read_csv(f_stitch,na_filter=False)
df_detect=pd.read_csv(f_detect,na_filter=False)

#Intersection and union of files
#intersection
df_merged=pd.merge(df_stitch,df_detect,how="inner",on=['time','region','season','var','minlat','minlon','maxlat','maxlon','centlat','centlon'])
#union
df_allmerge=pd.merge(df_stitch,df_detect,how='outer',on=['time','region','season','var','minlat','minlon','maxlat','maxlon','centlat','centlon'])
df_hasmerged=df_allmerge[df_allmerge['bnum_y'].isna()]
df_issplit=df_allmerge[df_allmerge['bnum_x'].isna()]

df_merged['bnum2'] = df_merged['bnum_x']
df_allmerge['bnum2'] = df_allmerge['bnum_x']
df_hasmerged['bnum2'] = df_hasmerged['bnum_x']
df_issplit['bnum2'] = df_issplit['bnum_y']   

df_total=pd.DataFrame()
for c in df_merged.columns.values:
    if not "_y" in c:
        if "_x" in c:
            newname=c.replace("_x","")
            df_total[newname]=df_merged[c]
        else:
            df_total[c]=df_merged[c]
df_total['merged']=['NO']*len(df_total)        
for t in sorted(list(set(df_hasmerged['time']))):
    df_checkmerge=df_hasmerged[df_hasmerged['time'].str.match(t)]
    df_splits=df_issplit[df_issplit['time'].str.match(t)]
    if (len(df_splits)>0):
        for n in range(0,len(df_checkmerge)):
            minlat_merged=df_checkmerge['minlat'].iloc[n]
            maxlat_merged=df_checkmerge['maxlat'].iloc[n]
            minlon_merged=df_checkmerge['minlon'].iloc[n]
            maxlon_merged=df_checkmerge['maxlon'].iloc[n]
            PER_BOUND=False
            if (minlon_merged>maxlon_merged):
                PER_BOUND=True
                minlon_merged=lon_convert_360_180(minlon_merged)
                maxlon_merged=lon_convert_360_180(maxlon_merged)
                if (minlon_merged>maxlon_merged):
                    maxlon_merged+=360.
            for m in range(0,len(df_splits)):
                minlat_split=df_splits['minlat'].iloc[m]
                maxlat_split=df_splits['maxlat'].iloc[m]
                minlon_split=df_splits['minlon'].iloc[m]
                maxlon_split=df_splits['maxlon'].iloc[m]


                minlon_split_180=lon_convert_360_180(minlon_split)
                maxlon_split_180=lon_convert_360_180(maxlon_split)

                #Check original boundaries
                is_inside=False
                if ((minlat_split>=minlat_merged) &
                    (maxlat_split<=maxlat_merged)&
                    (minlon_split>=minlon_merged)&
                    (maxlon_split<=maxlon_merged)):
                    is_inside=True
                else:
                    if ((minlat_split>=minlat_merged) &
                    (maxlat_split<=maxlat_merged)&
                    (minlon_split_180>=minlon_merged)&
                    (maxlon_split_180<=maxlon_merged)):
                        is_inside=True
                if (is_inside==True):
                    merge_dict={}
                    for c in df_splits.columns.values:
                        if not "_x" in c:
                            if "_y" in c:
                                if ((c!="bnum_y") & (c!="fname_y")):
                                    newname=c.replace("_y","")
                                    merge_dict[newname]=df_splits[c].iloc[m]
                            else:
                                merge_dict[c] = df_splits[c].iloc[m]
                    merge_dict['bnum']=df_checkmerge['bnum_x'].iloc[n]
                    merge_dict['merged']='YES'
                    merge_dict['fname']=df_checkmerge['fname_x'].iloc[n]
                    df_total=df_total.append(merge_dict,ignore_index=True)
df_total=df_total.sort_values(['var','region','time'])
#print("Length of stitchblobs for {:} {:} was {:} and length of detectblobs was {:}; length of merged dataframe is {:}".format(d,r,len(df_stitch),len(df_detect),len(df_total)))

df_total.to_csv(fname_out,index=False,na_rep="_")