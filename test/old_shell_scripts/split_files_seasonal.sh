#!/bin/bash

YEAR=""
SEASON=""
DATA=""
DIR=""
for i in "$@"; do
  case $i in 
    --year=*)
    YEAR="${i#*=}"
    ;;
    --season=*)
    SEASON="${i#*=}"
    ;;
    --data=*)
    DATA="${i#*=}"
    ;;
    --dir=*)
    DIR="${i#*=}"
    ;;
    *)
    
    ;;
  esac
done

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
else
  echo "Error: don't recognize season (DJF, MAM, JJA, or SON)."
  exit
fi



MINIT=$((MONTH-1))
YF=$(printf "%04d" $YINIT)
M=$(printf "%02d" $MONTH)

#Check that file hasn't already been split
FILE_CHECK=$(ls $DIR/$DATA/data/*$YF-$M-01*devs.nc 2>/dev/null | wc -l)

if [[ $FILE_CHECK -gt 0 ]]; then
  echo "File has already been split. Exiting."
  exit
fi


MF=$(printf "%02d" $MINIT)
SPLIT_FILE=$(ls $DIR/$DATA/data/*$YF-$MF*devs.nc | tail -1)

~/tempestextremes/bin/split_file --in $SPLIT_FILE --rename --vars DIPV,ADIPV,INT_ADIPV
