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

avgname=$(ls $DATA_DIR/*avg.nc)

#Create batchfiles for jobs
#AVG script
echo "#PBS -q serial">devbatch
echo "#PBS -l walltime=3:00:00">>devbatch
echo "#PBS -l vmem=10GB">>devbatch
echo " ">>devbatch
echo "cd $WKDIR">>devbatch
echo "./blockingDev --inlist integ_list.txt --avg $avgname">>devbatch
#qsub.serial devbatch

