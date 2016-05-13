#!/bin/bash

#Automate the calculation of the blobs files + stats

SEASONS=("MAM" "JJA" "SON" "DJF")
DATA=("climo" "2xCO2" "SSTplus2")
mstart=(3 6 9 12)

for d in ${DATA[@]}; do
  cd $SCRATCH/$d/data
  bdir="$SCRATCH/$d/blobs"
  if [ ! -e $bdir ]; then
    mkdir -p $bdir
  fi
  for y in {2..23}; do
    i=0
    for s in ${SEASONS[@]}; do
      yf=$(printf "%04d" $y)
      m=${mstart[i]}
      mf=$(printf "%02d" $m)

      if [ $m -eq 12 ]; then
        ls *$yf-$mf*devs.nc > bloblist
        yf=$(printf "%04d" $((y+1)))
        ls *$yf-0[12]*devs.nc >> bloblist
      else
        mstring=$(printf "*$yf-{%02d,%02d,%02d}*devs.nc" $m $((m+1)) $((m+2)))
        echo "ls $mstring" | sh > bloblist
      fi
      #run StitchBlobs
      echo "$s""_""$yf"
      blobsname="$bdir/$s""_""$yf""_blobs_""$d.nc"
      statsname="$bdir/$s""_""$yf""_stats_""$d.txt"
      ~/tempestextremes/bin/StitchBlobs --inlist bloblist --out $blobsname --var INT_ADIPV --outvar PV_BLOB --minsize 5 --mintime 40
      ~/tempestextremes/bin/BlobStats --infile $blobsname --outfile $statsname --invar PV_BLOB --out minlat,maxlat,minlon,maxlon,centlat,centlon,area


      i=$((i+1))
    done
  done
done

