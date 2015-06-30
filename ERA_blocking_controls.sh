#!/bin/bash

#This script will automate the PV blocking calculations for the input data list

export ystart=$1
export yend=$2
export yint=$3
export DATA_DIR="/media/mariellep/my_hd/ERA_data"
export AVG_TXT="ERA_avg_list.txt"
export DEV_TXT="ERA_dev_list.txt"

if [ $# -lt 3 ]; then
  echo "Forgot to provide year numbers!"
  exit
fi

if [ ! -d $DATA_DIR ]; then
  echo "Error! Check your file path."
  exit
fi

cd $DATA_DIR

if [ -e $AVG_TXT ]; then
  rm $AVG_TXT
fi

if [ -e $DEV_TXT ]; then
  rm $DEV_TXT
fi

#Search for files over specified range
ys=$((ystart))
ye=$((yend))
for ((y=$ys; y<=$ye; y++)); do
  for m in $(seq -f "%02g" 1 12); do
    infile="ERA_"$y"/ERA_"$y"_"$m"_vars.nc"
    outfile="ERA_"$y"/ERA_"$y"_"$m"_vars_integ.nc"
    ~/tempestextremes/blockingPV --in $infile --out $outfile --hpa
    echo $outfile >> $AVG_TXT
  done
done

avg_outfile="ERA_avg/ERA_"$ys"_"$ye"_dailyavg.nc"
~/tempestextremes/blockingAvg --inlist $AVG_TXT --out $avg_outfile
~/tempestextremes/blockingDevs --inlist $AVG_TXT --avg $avg_outfile

yi=$((yint))
for m in $(seq -f "%02g" 1 12); do
  devfile="ERA_"$yi"/ERA_"$yi"_"$m"_vars_integ_devs.nc"
  echo $devfile >> $DEV_TXT
done

blobsfile="ERA_"$yi"_blobs.nc"
~/tempestextremes/StitchBlobs --inlist $DEV_TXT --out $blobsfile --var INT_ADIPV --minsize 5 --mintime 12

~/tempestextremes/blockingDensity --in $blobsfile --var INT_ADIPVtag --out "ERA_"$yi"_density.nc"
