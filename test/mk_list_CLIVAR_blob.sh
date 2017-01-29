#!/bin/bash

#Automate the calculation of the blobs files + stats

SEASONS=("MAM" "JJA" "SON" "DJF")
DATA=("climo" "2xCO2" "SSTplus2")
mstart=(3 6 9 12)

#Addition of regional parameters
SECTOR=("NA" "NC" "NP" "SA" "SI" "SP")
LEFT_BOUND=(250 30 130 290 20 120)
RIGHT_BOUND=(50 150 270 40 140 310)
MIN_LAT=(25 25 25 -75 -75 -75)
MAX_LAT=(75 75 75 -25 -25 -25)

DAT_TYPE="PV"
for d in ${DATA[@]}; do
  cd $SCRATCH/$d/blobs

    for s in ${SEASONS[@]}; do
      yf=$(printf "%04d" $y)
      for a in ${SECTOR[@]}; do
       if [ "$DAT_TYPE" == "PV" ]; then        
         mstr=$SCRATCH/$d/blobs/$s"_*_"$a"_dens_"$d".nc"
         listname="$SCRATCH/$d/blobs/"$d"_"$s"_"$a"_all_blobs_bloblist"
         densname="$SCRATCH/$d/blobs/"$d"_"$s"_"$a"_avg_dens_time.nc"
       else
         mstr=$SCRATCH/$d/blobs/$s"_*_"$a"_zdens_"$d".nc"
         listname="$SCRATCH/$d/blobs/"$d"_"$s"_"$a"_all_zblobs_bloblist"
         densname="$SCRATCH/$d/blobs/"$d"_"$s"_"$a"_avg_zdens_time.nc"
       fi
         nces $mstr $densname
      done
    done
done

