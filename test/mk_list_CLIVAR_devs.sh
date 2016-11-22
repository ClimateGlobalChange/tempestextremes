#!/bin/bash

#Automate the calculation of the blobs files + stats

SEASONS=("MAM" "JJA" "SON" "DJF")
DATA=("climo" "2xCO2" "SSTplus2" "SSTplus2_2xCO2")
mstart=(3 6 9 12)

#Addition of regional parameters
SECTOR=("NA" "NC" "NP" "SA" "SI" "SP")

for d in ${DATA[@]}; do
  for s in ${SEASONS[@]}; do
    listname="$SCRATCH/$d/data/"$d"_"$s"_all_devs_bloblist"
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
      listname="$SCRATCH/$d/data/"$d"_"$s"_all_devs_bloblist"
      yf=$(printf "%04d" $y)
      m=${mstart[i]}
      mf=$(printf "%02d" $m)

      if [ $m -eq 12 ]; then
        ls *$yf-$mf*integ_devs.nc > bloblist
        yf=$(printf "%04d" $((y+1)))
        ls *$yf-0[12]*integ_devs.nc >> bloblist
      else
        mstring=$(printf "*$yf-{%02d,%02d,%02d}*integ_devs.nc" $m $((m+1)) $((m+2)))
        echo "ls $mstring" | sh > bloblist
      fi
      i=$((i+1))
      cat bloblist >> $listname
    done
  done
done

