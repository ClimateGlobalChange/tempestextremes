#!/bin/bash

#YEAR=""
#SEASON=""
#DATA=""
#DIR=""
#for i in "$@"; do
#  case $i in 
#    --year=*)
#    YEAR="${i#*=}"
#    ;;
#    --season=*)
#    SEASON="${i#*=}"
#    ;;
#    --data=*)
#    DATA="${i#*=}"
#    ;;
#    --dir=*)
#    DIR="${i#*=}"
#    ;;
#    *)
#    
#    ;;
#  esac
#done
DATAS=("climo" "2xCO2" "SSTplus2")
for d in ${DATAS[@]}; do
DIR=$SCRATCH/$d
SEASONS=("DJF" "MAM" "JJA" "SON")
for s in ${SEASONS[@]}; do
  SEASON=$s
for YEAR in {2..24}; do
  YINIT=$((YEAR))
  MONTH=0
  if [[ "$SEASON" == "DJF" ]]; then
    YINIT=$((YEAR-1))
    MONTH=12  
  elif [[ "$SEASON" == "MAM" ]]; then
    MONTH=3
  elif [[ "$SEASON" == "JJA" ]]; then
    MONTH=6
  elif [[ "$SEASON" == "SON" ]]; then
    MONTH=9
  fi



  MINIT=$((MONTH-1))
  YF=$(printf "%04d" $YINIT)
  M=$(printf "%02d" $MONTH)

#Check that file hasn't already been split
  FILE_CHECK=$(ls $DIR/data/*$YF-$M-01*integ.nc 2>/dev/null | wc -l)

  if [[ $FILE_CHECK -gt 0 ]]; then
   echo "File has already been split."
  else


    MF=$(printf "%02d" $MINIT)
    SPLIT_FILE=$(ls $DIR/data/*$YF-$MF*integ.nc | tail -1)
    if [[ "$SPLIT_FILE" == "" ]]; then
      echo "File doesn't exist."
    else
      echo $SPLIT_FILE
      ~/tempestextremes/bin/split_file --in $SPLIT_FILE --rename --vars IPV,AVGT,AVGU,AVGV
    fi
  fi
done
done
done
