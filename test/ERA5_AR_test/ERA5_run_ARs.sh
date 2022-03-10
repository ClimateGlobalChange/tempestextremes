#!/bin/bash

TEMPESTEXTREMESDIR=/global/homes/p/paullric/tempestextremes/bin

mkdir -p ERA5_AR

rm -rf ERA5_IVT_files.txt
rm -rf ERA5_AR_files.txt
rm -rf ERA5_AR_SB_files.txt
for f in /global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.vinteg/1979*
do
        yearmonth=$(basename $f)
        for g in $f/*162_071_viwve*
        do
                h=$(basename $g)
                suffix=${h: -33}
                prefix=${h:0:18}
                echo "/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.vinteg/${yearmonth}/${prefix}162_071_viwve${suffix};/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.vinteg/${yearmonth}/${prefix}162_072_viwvn${suffix}" >> ERA5_IVT_files.txt
                echo "ERA5_AR/e5.oper.an.vinteg.ar_tag${suffix}" >> ERA5_AR_files.txt
                echo "ERA5_AR/e5.oper.an.vinteg.ar_tag_stitch${suffix}" >> ERA5_AR_SB_files.txt
        done
done

srun -n 32 $TEMPESTEXTREMESDIR/DetectBlobs --in_data_list ERA5_IVT_files.txt --out_list ERA5_AR_files.txt --thresholdcmd "_LAPLACIAN{8,10.0}(_VECMAG(VIWVE,VIWVN)),<=,-30000,0" --minabslat 20 --geofiltercmd "area,>,850000km2" --latname "latitude" --lonname "longitude" --timefilter "6hr" --tagvar "AR_binary_tag"

