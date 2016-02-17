#!/bin/bash

#This script will automate the PV blocking calculations for the input data list

export DATA_DIR=$GSCRATCH/CO2_2
export WKDIR=$HOME/tempestextremes
export NUM_NODES=8
if [ ! -d $DATA_DIR ]; then
  echo "Error! Check your file path."
  exit
fi

#Create batchfiles for jobs
#PV script
for (( i=1; i<=$NUM_NODES; i++ )); do
  fname="pvbatch$i"
  echo "#PBS -q serial">$fname
  echo "#PBS -l walltime=3:00:00">>$fname
  echo "#PBS -l vmem=10GB">>$fname
  echo " ">>$fname
  echo "cd $WKDIR">>$fname
done

#Check if files have already been interpolated
ipl_files=$(ls $DATA_DIR/h4/*_ipl.nc 2> /dev/null | wc -l)
if [ "$ipl_files" != "0" ]; then
  #Divide up jobs among nodes
  div=$(($ipl_files/$NUM_NODES))
  rem=$(($ipl_files%$NUM_NODES))
  ls $DATA_DIR/h4/*_ipl.nc > ipllist.txt
else
  ls $DATA_DIR/h2/*00000.nc > h2list.txt
  ls $DATA_DIR/h4/*00000.nc > h4list.txt
  paste h2list.txt h4list.txt > hlist.txt
  num_files=$(ls $DATA_DIR/h2/*.nc 2> /dev/null | wc -l)
  div=$(($num_files/$NUM_NODES))
  rem=$(($num_files%$NUM_NODES))
fi  
#write files and submit to queue
for (( i=1; i<=$NUM_NODES; i++ )); do
  #start and end lines of file
  startline=$((1+($i-1)*$div))
  endline=$(($i*$div))
  #batch name
  fname="pvbatch$i"
  if [ $i -eq $NUM_NODES ]; then
    endline=$((endline+$rem))
  fi
  if [ "$ipl_files" != "0" ]; then
    for (( x=$startline; x<=$endline; x++ )); do
      integ_line=$(sed -n "${x}p" "ipllist.txt")
      echo "./blockingPV --in $integ_line">>$fname
    done
  else
    for (( x=$startline; x<=$endline; x++ )); do
      readline=$(sed -n "${x}p" "hlist.txt")
      read h2 h4 <<<"$readline"
      echo "./blockingPV --in2d $h2 --in $h4 --ipl">>$fname
    done
  fi
   qsub.serial $fname
done

 

