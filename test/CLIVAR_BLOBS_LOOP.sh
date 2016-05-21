#!/bin/bash

#Automate the calculation of the blobs files + stats

SEASONS=("MAM" "JJA" "SON" "DJF")
DATA=("climo" "2xCO2" "SSTplus2" "SSTplus2_2xCO2")
mstart=(3 6 9 12)
BLOBS="FALSE"

if [ "$BLOBS" == "TRUE" ]; then
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
      densname="$bdir/$s""_""$yf""_dens_""$d.nc"
      ~/tempestextremes/bin/StitchBlobs --inlist bloblist --out $blobsname --var INT_ADIPV --outvar PV_BLOB --minsize 5 --mintime 40
      ~/tempestextremes/bin/BlobStats --infile $blobsname --outfile $statsname --invar PV_BLOB --out minlat,maxlat,minlon,maxlon,centlat,centlon,area
      ~/tempestextremes/bin/DensityCalculations --in $blobsname --var PV_BLOB --out $densname

      i=$((i+1))
    done
  done
done

#removing the files that contain missing dates!
#climo dataset
rm $SCRATCH/climo/blobs/DJF_0004_*
rm $SCRATCH/climo/blobs/DJF_0010_*
rm $SCRATCH/climo/blobs/JJA_0012_*
rm $SCRATCH/climo/blobs/SON_0020_*

#2xCO2 dataset
rm $SCRATCH/2xCO2/blobs/SON_0007_*

#SSTplus2 dataset
rm $SCRATCH/SSTplus2/blobs/SON_0016_*

#SSTplus2_2xCO2 dataset
rm $SCRATCH/SSTplus2_2xCO2/blobs/DJF_0007_*
rm $SCRATCH/SSTplus2_2xCO2/blobs/DJF_0016_*

fi


#Now run average density calculations
for d in ${DATA[@]}; do
  cd $SCRATCH/$d/blobs
  for s in ${SEASONS[@]}; do
    lsname=$s"_blobs"
    outname=$s"_avg_dens_"$d".nc"
    ls $s*blobs*.nc > $lsname
    numfiles=$(cat $lsname | wc -l)
    #~/tempestextremes/bin/DensityCalculations --std --inlist $lsname --var PV_BLOB --out $outname
    python ~/tempestextremes/test/plot_density.py $outname "$d $numfiles yr $s" "avg" "PV blocking"
  done
done
