#!/bin/bash

SEASONS=("MAM" "JJA" "SON" "DJF")
mstart=(3 6 9 12)
ystart=1980
ycalc=2005
yend=2000
DDIR=/Volumes/ExFAT_drive/ERA_files
BDIR=$DDIR/ERA_blobs
VAR="GH"

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





if [ ! -e $BDIR ]; then
  mkdir -p $BDIR
fi

for ((y=ystart; y<=ycalc; y++)); do
  i=0
  SUBDIR=$DDIR/ERA_$y
  echo "Entering subdirectory $SUBDIR"
  if [ ! -e $SUBDIR ]; then
    echo "Missing directory. Check connection."
    exit
  fi
  for s in ${SEASONS[@]}; do

    if [ "$s" == "DJF" ]; then
      echo "DJF!"
      mstring="$SUBDIR/ERA*_12_vars_$SUFF"
      echo "ls $mstring" | sh > bloblist

      yn=$((y+1)) 
      SUBDIR2=$DDIR/ERA_$yn

      mstring="$SUBDIR2/ERA*_0[12]_vars_$SUFF"
      echo "ls $mstring" | sh >>bloblist
     # if [ "$VAR" == "IPV" ]; then
     #   ls $SUBDIR/*_12_vars_integ_devs.nc > bloblist
     # fi
     # else if [ "$VAR" == "GH" ]; then
     #   ls $SUBDIR/*_12_vars_GH_devs.nc > bloblist
     # fi

     # if [ "$VAR" == "IPV" ]; then
     #   ls $SUBDIR2/*_0[12]_vars_integ_devs.nc >> bloblist
     # fi
     # else if [ "$VAR" == "GH" ]; then
     #   ls $SUBDIR2/*_0[12]_vars_GH_devs.nc >> bloblist
     # fi

    else
     m=${mstart[i]}
     # m1=${mstart[i]}
     # m2=$((m1+2))
     # echo "m1 is $m1 and m2 is $m2"
     # echo "ls $SUBDIR/*_*{$m1..$m2}_*devs.nc" | sh > bloblist
     
     mstring=$(printf "$SUBDIR/ERA*_{%02d,%02d,%02d}_vars_$SUFF" $m $((m+1)) $((m+2)))
     echo $mstring
     echo "ls $mstring" | sh > bloblist
    fi
    cat bloblist
    
    blobsname="$BDIR/ERA_"$y"_"$s"_"$BLOB_SUFF
    statsname="$BDIR/ERA_"$y"_"$s"_"$STAT_SUFF
    densname="$BDIR/ERA_"$y"_"$s"_"$DENS_SUFF
    vdensname="$BDIR/ERA_"$y"_"$s"_var_"$DENS_SUFF

    ~/tempestextremes/bin/StitchBlobs --inlist bloblist --out $blobsname --var $INVAR --outvar $BLOB_VAR --mintime 20
    ~/tempestextremes/bin/BlobStats --infile $blobsname --outfile $statsname --invar $BLOB_VAR --out minlat,maxlat,minlon,maxlon,centlat,centlon,area
    ~/tempestextremes/bin/DensityCalculations --in $blobsname --var $BLOB_VAR --out $densname
    ~/tempestextremes/bin/DensityCalculations --inlist bloblist --var $INVAR --out $vdensname
    i=$((i+1))
  done
done

cd $BDIR
for s in ${SEASONS[@]}; do
  lsname="ERA_"$s"_blobs"
  if [ -e $lsname ]; then
    rm $lsname
  fi
  outname="$BDIR/ERA_"$s"_avg_"$DENS_SUFF
  for ((y=ystart; y<=yend; y++)); do
    echo "ls ERA_"$y"_"$s"_"$BLOB_SUFF" >> $lsname" | sh
  done
  cat $lsname

  numfiles=$(cat $lsname | wc -l)
  n=$((numfiles))
  echo "There are $n files"
  ~/tempestextremes/bin/DensityCalculations --std --inlist $lsname --var $BLOB_VAR --out $outname
  python ~/tempestextremes/test/plot_density.py $outname "ERA $n yr $s" "avg" "$PLOT_TITLE"
done

cp /Volumes/ExFAT_drive/ERA_files/ERA_blobs/ERA*plot.png ~/figs/
