#!/bin/bash

ystart=$1
yend=$2
DATA_DIR="/media/mariellep/my_hd/ERA_data"
SEASONS=( "DJF" "MAM" "JJA" "SON" )
NDAYS=5


if [ $# -lt 2 ]; then
  echo "Need to provide start and end year numbers!"
  exit
fi

if [ ! -d $DATA_DIR ]; then 
  echo "Error! Check your file path."
  exit
fi

cd $DATA_DIR


if [ -e $BLOB_INPUT ]; then
  rm $BLOB_INPUT
fi

ys=$((ystart))
ye=$((yend))
nsteps=$((NDAYS*4))

#file for Dec of prev year
yprev=$((ys-1))
dfile="ERA_"$yprev"/ERA_"$yprev"_12_vars_GHblock.nc"

#Run TM method, StitchBlobs, Density for year and DJF
for ((y=$ys; y<=$ye; y++)); do
  BLOBLIST="ERA_blob_"$y"_GH_list.txt"
  if [ -e $BLOBLIST ]; then
    rm $BLOBLIST
  fi
  for m in $(seq -f "%02g" 1 12); do
    infile="ERA_"$y"/ERA_"$y"_"$m"_vars.nc"
    outfile="ERA_"$y"/ERA_"$y"_"$m"_vars_GHblock.nc"
    if [ ! -e $outfile ]; then
      ~/tempestextremes/blockingGH --in $infile --out $outfile --hpa --is4d
    fi
    echo $outfile >> $BLOBLIST
  done
  blobsfile="ERA_"$y"_"$NDAYS"day_GHblobs.nc"
  ~/tempestextremes/StitchBlobs --inlist $BLOBLIST --out $blobsfile --var GHGrad --minsize 5 --mintime $nsteps
  densfile="ERA_"$y"_GHdensity.nc"
  ~/tempestextremes/blockingDensity --in $blobsfile --var GHGradtag --out $densfile
  python ~/tempestextremes/plot_density.py $densfile year $y "GH gradient"

done

