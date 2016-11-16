#!/bin/bash

SEASONS=("MAM" "JJA" "SON" "DJF")
mstart=(3 6 9 12)
ystart=1980
ycalc=2005
yend=2000
DDIR=/Volumes/ExFAT_drive/ERA_files
BDIR=$DDIR/ERA_blobs
VAR="IPV"

#Addition of regional parameters
SECTOR=("NA" "NC" "NP" "SA" "SI" "SP")
LEFT_BOUND=(250 30 130 290 20 120)
RIGHT_BOUND=(50 150 270 40 140 310)
MIN_LAT=(25 25 25 -75 -75 -75)
MAX_LAT=(75 75 75 -25 -25 -25)

SUFF=""  
BLOB_SUFF=""  
STAT_SUFF=""
DENS_SUFF=""
INVAR=""
BLOB_VAR=""
PLOT_TITLE=""
if [ "$VAR" == "IPV" ]; then
  SUFF="integ_devs.nc"
  BLOB_SUFF="blobs.nc"
  STAT_SUFF="stats.txt"
  DENS_SUFF="dens.nc"
  INVAR="INT_ADIPV"
  BLOB_VAR="PV_BLOB"
  PLOT_TITLE="PV blocking"
elif [ "$VAR" == "GH" ]; then
  SUFF="GH_devs.nc"
  BLOB_SUFF="GHblobs.nc"
  STAT_SUFF="GHstats.txt"
  DENS_SUFF="GHdens.nc"
  INVAR="INT_ADGH"
  BLOB_VAR="GH_BLOB"
  PLOT_TITLE="GH blocking"
fi


SUFF="integ.nc"


for s in ${SEASONS[@]}; do
  listname=$DDIR"/ERA_"$s"_all_bloblist"
  if [[ -e $listname ]]; then
    rm $listname
  fi
done

for ((y=ystart; y<=ycalc; y++)); do
  i=0
  SUBDIR=$DDIR/ERA_$y
  echo "Entering subdirectory $SUBDIR"
  if [ ! -e $SUBDIR ]; then
    echo "Missing directory. Check connection."
    exit
  fi
  for s in ${SEASONS[@]}; do
    listname=$DDIR"/ERA_"$s"_all_bloblist"
    if [ "$s" == "DJF" ]; then
      mstring="$SUBDIR/ERA*_12_vars_$SUFF"
      echo "ls $mstring" | sh > bloblist

      yn=$((y+1)) 
      SUBDIR2=$DDIR/ERA_$yn

      mstring="$SUBDIR2/ERA*_0[12]_vars_$SUFF"
      echo "ls $mstring" | sh >>bloblist

    else
     m=${mstart[i]}
     
     mstring=$(printf "$SUBDIR/ERA*_{%02d,%02d,%02d}_vars_$SUFF" $m $((m+1)) $((m+2)))
     echo $mstring
     echo "ls $mstring" | sh > bloblist
    fi
    cat bloblist >> $listname


    i=$((i+1))
  done
done

