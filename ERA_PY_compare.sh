#!/bin/bash

#This script will automate the PV blocking calculations for the input data list

export DATA_DIR="/media/mariellep/my_hd/ERA_data/PV_anomaly_calcs/PY_ERA_CALC/"
export BLOB_INPUT="ERA_py_blob_list.txt"
export SEASONS=( "DJF" "MAM" "JJA" "SON" )
export NDAYS=5


if [ ! -d $DATA_DIR ]; then
  echo "Error! Check your file path."
  exit
fi

cd $DATA_DIR


if [ -e $BLOB_INPUT ]; then
  rm $BLOB_INPUT
fi



c=0
#Make list to input for StitchBlobs
for m in $(seq -f "%02g" 1 12); do
  devfile="ERA_devs/ERA_2013_"$m"_1981to2005_intdevs.nc"
  if [ ! -e $devfile ]; then
    c=$((c+1))
  fi
  echo $devfile >> $BLOB_INPUT
done

echo $c " files missing"

nsteps=$((NDAYS*4))
echo $nsteps

yi=2013

#Run StitchBlobs for whole year
blobsfile="ERA_"$yi"_"$NDAYS"day_blobs.nc"
~/tempestextremes/StitchBlobs --inlist $BLOB_INPUT --out $blobsfile --var INT_ADIPV --minsize 5 --mintime $nsteps

#Calculate blocking density for year
densfile="ERA_"$yi"_density.nc"
~/tempestextremes/blockingDensity --in $blobsfile --var INT_ADIPVtag --out $densfile

#Make plot of density for year
python ~/tempestextremes/plot_density.py $densfile year $yi 

#Run StitchBlobs for seasonal (DJF,MAM,JJA,SON)
#DJF
#check that December (previous year!) file exists
dfile="ERA_devs/ERA_2012_12_1981to2005_intdevs.nc"
if [ ! -e $dfile ]; then
  echo "Error: deviations not calculated for December of "$yprev
  exit
else
  #array of season list file names
  x=0
  declare -a SFILES
  for s in ${SEASONS[@]}; do
    echo $x
    SFILES[x]="ERA_"$s"_blobs_list.txt"
    x=$((x+1))
  done
  #initialize lists
  echo $dfile > ${SFILES[0]}
  echo "ERA_devs/ERA_"$yi"_03_1981to2005_intdevs.nc" > ${SFILES[1]}
  echo "ERA_devs/ERA_"$yi"_06_1981to2005_intdevs.nc" > ${SFILES[2]}
  echo "ERA_devs/ERA_"$yi"_09_1981to2005_intdevs.nc" > ${SFILES[3]}
  #other files added to lists
  for n in $(seq 1 2); do
    for x in $(seq 0 3); do
      n1=$((3*$x+$n))
      n2=$(printf "%02d" $n1)
      echo "ERA_devs/ERA_"$yi"_"$n2"_1981to2005_intdevs.nc" >> ${SFILES[$x]}
    done
  done
fi

# Run StitchBlobs and density calcs for seasonal
for x in $(seq 0 3); do 
  blobsfile="ERA_"$yi"_"${SEASONS[$x]}"_"$NDAYS"day_blobs.nc"
  densfile="ERA_"$yi"_"${SEASONS[$x]}"_density.nc"

  ~/tempestextremes/StitchBlobs --inlist ${SFILES[$x]} --out $blobsfile --var INT_ADIPV --minsize 5 --mintime $nsteps
  ~/tempestextremes/blockingDensity --in $blobsfile --var INT_ADIPVtag --out $densfile
  python ~/tempestextremes/plot_density.py $densfile ${SEASONS[$x]} $yi 
done

