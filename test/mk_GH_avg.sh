DATAS=("climo" "SSTplus2" "2xCO2")
SEASONS=("DJF" "MAM" "JJA" "SON")


for d in ${DATAS[@]}; do
  cd $SCRATCH/$d/data
  for s in ${SEASONS[@]}; do
    inlist=$d"_"$s"_all_blobzlist"
    devlist=$d"_"$s"_all_devs_blobzlist"
    outname=$d"_"$s"_avg_Z_time.nc"
    outdev=$d"_"$s"_avg_Z_devs_time.nc"
    ~/tempestextremes/bin/AvgVar --inlist $inlist --varlist Z3 --out $outname
    ~/tempestextremes/bin/AvgVar --inlist $devlist --varlist DGH,ADGH,INT_ADGH --out $outdev
  done
done
