#!/bin/bash

#This script will automate the PV blocking calculations for the input data list

export ystart=$1
export yend=$2
export yint=$3
export DATA_DIR="/media/mariellep/my_hd/ERA_data"
export AVG_INPUT="ERA_avg_list.txt"
export DEV_INPUT="ERA_dev_list.txt"
export BLOB_INPUT="ERA_blob_list.txt"

if [ $# -lt 3 ]; then
  echo "Forgot to provide year numbers!"
  exit
fi

if [ ! -d $DATA_DIR ]; then
  echo "Error! Check your file path."
  exit
fi

cd $DATA_DIR

if [ -e $AVG_INPUT ]; then
  rm $AVG_INPUT
fi

if [ -e $DEV_INPUT ]; then
  rm $DEV_INPUT
fi

if [ -e $BLOB_INPUT ]; then
  rm $BLOB_INPUT
fi

#Search for files over specified range
ys=$((ystart))
ye=$((yend))
for ((y=$ys; y<=$ye; y++)); do
  for m in $(seq -f "%02g" 1 12); do
    infile="ERA_"$y"/ERA_"$y"_"$m"_vars.nc"
    outfile="ERA_"$y"/ERA_"$y"_"$m"_vars_integ.nc"
    #run integration code if file doesn't already exist
    if [ !-e $outfile ]; then
      ~/tempestextremes/blockingPV --in $infile --out $outfile --hpa
    fi
    echo $outfile >> $AVG_INPUT
  done
done

#List of files for year that PV* is calculated
yi=$((yint))
for m in $(seq -f "%02g" 1 12); do
  infile="ERA_"$yi"/ERA_"$yi"_"$m"_vars_integ.nc"
  echo $infile >> $DEV_INPUT
done

#Calculate block average
avg_outfile="ERA_avg/ERA_"$ys"_"$ye"_dailyavg.nc"
if [ !-e $avg_outfile ]; then
  ~/tempestextremes/blockingAvg --inlist $AVG_INPUT --out $avg_outfile
fi

#calculate deviations
~/tempestextremes/blockingDevs --inlist $DEV_INPUT --avg $avg_outfile

#Make list to input for StitchBlobs
for m in $(seq -f "%02g" 1 12); do
  devfile="ERA_"$yi"/ERA_"$yi"_"$m"_vars_integ_devs.nc"
  echo $devfile >> $BLOB_INPUT
done

#Run StitchBlobs
blobsfile="ERA_"$yi"_blobs.nc"
~/tempestextremes/StitchBlobs --inlist $BLOB_INPUT --out $blobsfile --var INT_ADIPV --minsize 5 --mintime 12

#Calculate blocking density
~/tempestextremes/blockingDensity --in $blobsfile --var INT_ADIPVtag --out "ERA_"$yi"_density.nc"
