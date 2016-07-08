#!/bin/bash

#This script will automate the PV blocking calculations for the input data list

export DATA_DIR=$GSCRATCH/CO2_2/h4
export WKDIR=$HOME/tempestextremes
export AVG_TXT="cam_avg_list.txt"

if [ ! -d $DATA_DIR ]; then
  echo "Error! Check your file path."
  exit
fi

if [ ! -e integ_list.txt ]; then
  ls $DATA_DIR/*integ.nc > integ_list.txt
fi

fname="$(cat integ_list.txt | head -n 1)"
avgname="$(echo $fname | sed s/0.*integ/avg/)"

#Create batchfiles for jobs
#AVG script
echo "#PBS -q serial">avgbatch
echo "#PBS -l walltime=3:00:00">>avgbatch
echo "#PBS -l vmem=10GB">>avgbatch
echo " ">>avgbatch
echo "cd $WKDIR">>avgbatch
echo "./blockingAvg --inlist integ_list.txt --out $avgname">>avgbatch
qsub.serial avgbatch

#~/tempestextremes/blockingAvg --inlist $AVG_TXT --out $avg_outfile
#~/tempestextremes/blockingDevs --inlist $AVG_TXT --avg $avg_outfile

#yi=$((yint))
#for m in $(seq -f "%02g" 1 12); do
#  devfile="ERA_"$yi"/ERA_"$yi"_"$m"_vars_integ_devs.nc"
#  echo $devfile >> $DEV_TXT
#done

#blobsfile="ERA_"$yi"_blobs.nc"
#~/tempestextremes/StitchBlobs --inlist $DEV_TXT --out $blobsfile --var INT_ADIPV --minsize 5 --mintime 12

#~/tempestextremes/blockingDensity --in $blobsfile --var INT_ADIPVtag --out "ERA_"$yi"_density.nc"
