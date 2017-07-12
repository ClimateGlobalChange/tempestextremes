#!/bin/bash

SEASONS=("MAM" "JJA" "SON" "DJF")
mstart=(3 6 9 12)
ystart=1980
ycalc=2003
yend=2000
DDIR=/Volumes/ExFAT_drive/ERA_files
BDIR=$DDIR/ERA_detect
#VAR="GHGrad"

VARVEC=("Z" "IPV")

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

for VAR in ${VARVEC[@]}; do

if [ "$VAR" == "IPV" ]; then
  SUFF="integ_devs_norm_const.nc"
  BLOB_SUFF="blobs_nostitch_const.nc"
  STAT_SUFF="stats_nostitch_const.txt"
  DENS_SUFF="dens_nostitch_const.nc"
  INVAR="INT_ADIPV"
  BLOB_VAR="PV_BLOB"
  PLOT_TITLE="PV blocking"
elif [ "$VAR" == "Z" ]; then
  SUFF="z500_devs_norm_const.nc"
  BLOB_SUFF="Zblobs_nostitch_const.nc"
  STAT_SUFF="Zstats_nostitch_const.txt"
  DENS_SUFF="Zdens_nostitch_const.nc"
  INVAR="INT_ADGH"
  BLOB_VAR="Z_BLOB"
  PLOT_TITLE="Z blocking"
elif [ "$VAR" == "GHGrad" ]; then
  SUFF="z500_GHG.nc"
  BLOB_SUFF="GHGblobs_nostitch.nc"
  STAT_SUFF="GHGstats_nostitch.txt"
  DENS_SUFF="GHGdens_nostitch.nc"
  INVAR="GHGrad"
  BLOB_VAR="GHG_BLOB"
  PLOT_TITLE="GHG blocking"
fi



if [ ! -e $BDIR ]; then
  mkdir -p $BDIR
fi

for ((y=1980; y<=2005; y++)); do
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
    for n in {0..5}; do
      secname=${SECTOR[n]}
#    blobsname="$BDIR/ERA_"$y"_"$s"_"$BLOB_SUFF
#    statsname="$BDIR/ERA_"$y"_"$s"_"$STAT_SUFF
#    densname="$BDIR/ERA_"$y"_"$s"_"$DENS_SUFF
#    vdensname="$BDIR/ERA_"$y"_"$s"_var_"$DENS_SUFF

      blobsname="$BDIR/ERA_"$y"_"$s"_""$secname""_"$BLOB_SUFF
      statsname="$BDIR/ERA_"$y"_"$s"_""$secname""_"$STAT_SUFF
      densname="$BDIR/ERA_"$y"_"$s"_""$secname""_"$DENS_SUFF
      vdensname="$BDIR/ERA_"$y"_"$s"_var_""$secname""_"$DENS_SUFF

      ~/tempestextremes/bin/DetectBlobs --inlist bloblist --out $blobsname --var $INVAR --outvar $BLOB_VAR --minlat ${MIN_LAT[n]} --maxlat ${MAX_LAT[n]} --minlon ${LEFT_BOUND[n]} --maxlon ${RIGHT_BOUND[n]} --thresholdcmd "minarea,1000000000000"
      ~/tempestextremes/bin/BlobStats --infile $blobsname --outfile $statsname --invar $BLOB_VAR --out minlat,maxlat,minlon,maxlon,centlat,centlon,area
   #   ~/tempestextremes/bin/DensityCalculations --in $blobsname --var $BLOB_VAR --out $densname
   #   ~/tempestextremes/bin/DensityCalculations --inlist bloblist --var $INVAR --out $vdensname
      n=$((n+1))
    done
    i=$((i+1))
  done
done

#cd $BDIR
#for s in ${SEASONS[@]}; do
#  lsname="ERA_"$s"_blobs"
#  if [ -e $lsname ]; then
#    rm $lsname
#  fi
#  outname="$BDIR/ERA_"$s"_avg_"$DENS_SUFF
#  for ((y=ystart; y<=yend; y++)); do
#    echo "ls ERA_"$y"_"$s"_"$BLOB_SUFF" >> $lsname" | sh
#  done
#  cat $lsname

#  numfiles=$(cat $lsname | wc -l)
#  n=$((numfiles))
#  echo "There are $n files"
#  ~/tempestextremes/bin/DensityCalculations --std --inlist $lsname --var $BLOB_VAR --out $outname
#  python ~/tempestextremes/test/plot_density.py $outname "ERA $n yr $s" "avg" "$PLOT_TITLE"
#done
done
#cp /Volumes/ExFAT_drive/ERA_files/ERA_blobs/ERA*plot.png ~/figs/
