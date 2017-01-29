#!/bin/bash

#Automate the calculation of the blobs files + stats

SEASONS=("MAM" "JJA" "SON" "DJF")
DATA=("climo" "2xCO2" "SSTplus2" "SSTplus2_2xCO2")
mstart=(3 6 9 12)

#Addition of regional parameters
SECTOR=("NA" "NC" "NP" "SA" "SI" "SP")
LEFT_BOUND=(250 30 130 290 20 120)
RIGHT_BOUND=(50 150 270 40 140 310)
MIN_LAT=(25 25 25 -75 -75 -75)
MAX_LAT=(75 75 75 -25 -25 -25)

LIST_TYPE="GH"

for d in ${DATA[@]}; do
  for s in ${SEASONS[@]}; do
    if [ "$LIST_TYPE" == "GH" ]; then
      listname="$SCRATCH/$d/data/"$d"_"$s"_all_blobzlist"
    else
      listname="$SCRATCH/$d/data/"$d"_"$s"_all_devs_blobzlist"
    fi
    if [[ -e $listname ]]; then
      rm $listname
    fi
  done
done

for d in ${DATA[@]}; do
  cd $SCRATCH/$d/data

  for y in {2..23}; do
    i=0
    for s in ${SEASONS[@]}; do
      if [ "$LIST_TYPE" == "GH" ]; then      
        listname="$SCRATCH/$d/data/"$d"_"$s"_all_blobzlist"
        suffix="z500"
      else
        listname="$SCRATCH/$d/data/"$d"_"$s"_all_devs_blobzlist"
        suffix="z500_devs"
      fi
      yf=$(printf "%04d" $y)
      m=${mstart[i]}
      mf=$(printf "%02d" $m)
      if [ $m -eq 12 ]; then
        ls *$yf-$mf*$suffix.nc > blobzlist
        yf=$(printf "%04d" $((y+1)))
        ls *$yf-0[12]*$suffix.nc >> blobzlist
      else
        mstring=$(printf "*$yf-{%02d,%02d,%02d}*$suffix.nc" $m $((m+1)) $((m+2)))
        echo "ls $mstring" | sh > blobzlist
      fi
      i=$((i+1))
      cat blobzlist >> $listname
    done
  done
done

