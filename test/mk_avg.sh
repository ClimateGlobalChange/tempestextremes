DATAS=("climo" "SSTplus2" "2xCO2")
SEASONS=("DJF" "MAM" "JJA" "SON")


for d in ${DATAS[@]}; do
  cd $SCRATCH/$d/data
  for s in ${SEASONS[@]}; do
    inlist=$d"_"$s"_all_bloblist"
    devlist=$d"_"$s"_all_devs_bloblist"
    outname=$d"_"$s"_avg_time.nc"
    outdev=$d"_"$s"_avg_devs_time.nc"
#    ~/tempestextremes/bin/AvgVar --inlist $inlist --varlist IPV,AVGT,AVGU,AVGV --out $outname
    ~/tempestextremes/bin/AvgVar --inlist $devlist --varlist DIPV,ADIPV,INT_ADIPV --out $outdev
  done
done
