DATAS=("climo" "SSTplus2" "2xCO2")
SEASONS=("DJF" "MAM" "JJA" "SON")
SECTORS=("NA" "NC" "NP" "SA" "SI" "SP")

for d in ${DATAS[@]}; do
  cd $SCRATCH/$d/data
  for s in ${SEASONS[@]}; do
    inlist=$d"_"$s"_all_bloblist"
    devlist=$d"_"$s"_all_devs_bloblist"
    outname=$d"_"$s"_avg_time.nc"
    outdev=$d"_"$s"_avg_devs_time.nc"
    #~/tempestextremes/bin/AvgVar --inlist $inlist --varlist IPV,AVGT,AVGU,AVGV --out $outname
    #~/tempestextremes/bin/AvgVar --inlist $devlist --varlist DIPV,ADIPV,INT_ADIPV --out $outdev

  done
done

for d in ${DATAS[@]}; do
  cd $SCRATCH/$d/blobs
  for s in ${SEASONS[@]}; do
    for a in ${SECTORS[@]}; do
       inlist=$d"_"$s"_"$a"_blobs"
       outname=$d"_"$s"_"$a"_avg_blobs_time.nc"

       ~/tempestextremes/bin/AvgVar --inlist $inlist --varlist dens --out $outname
    done
  done
done
