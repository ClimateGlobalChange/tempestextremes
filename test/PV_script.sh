#!/bin/bash

# PV_script.sh
# Version 2.0 Mar 18, 2016 
# Author: Marielle Pinheiro

# This script automates PV calculations for input data
#  (specifically the cam5 runs)
# You can run the entire process (calculating PV, 
#  averages, deviations, and StitchBlobs) or only a
#  single job by specifying with the appropriate flags

#TO RUN THE SCRIPT:

#At the command line enter:
# ./PV_script.sh   WITH THE ADDITIONAL FLAGS:
# --ddir=[PATH WHERE DATA WILL BE STORED] ##MANDATORY
# --bdir=[PATH WHERE BINARIES ARE LOCATED] ##MANDATORY
# JOB FLAGS:
# --pv: Take h2 and h4 files and interpolate data to
#       pressure levels, then calculate potential vorticity
# --avg: find long term daily average using provided list of PV files
# --dev: find the deviations from the long term daily average
# FOR PV CALCULATIONS ONLY:
# --h2dir=[PATH WHERE H2 FILES LOCATED] ##MANDATORY
# --h4dir=[PATH WHERE H4 FILES LOCATED] ##MANDATORY
# --list2d=[LIST OF H2 FILES] ##OPTIONAL
# --list4d=[LIST OF H4 FILES] ##OPTIONAL
# --listpv=[LIST OF OUTPUT PV FILES] ##OPTIONAL
# FOR AVERAGE CALCULATIONS:
# --listpv=[LIST OF INPUT FILES] ##MANDATORY IF NOT RUNNING PV
# --varname=[NAME OF VARIABLE BEING AVERAGED] ##MANDATORY
# --avgname=[NAME OF AVERAGE VARIABLE] ##MANDATORY
# FOR DEVIATION CALCULATIONS:
# --listpv=[LIST OF INPUT FILES] ##MANDATORY IF NOT RUNNING PV
# --avgfile=[NAME OF AVERAGED FILE] ##MANDATORY IF NOT RUNNING AVG
#PLANNED UPDATES:
# Working on building in support for missing files
# Add to PV: add option for if interpolated files are already created
# Add StitchBlobs component to chained jobs
# Add GH anomaly 

DATA_DIR=""
BINDIR=""
H2DIR=""
H4DIR=""
LIST_4D=""
LIST_2D=""
LIST_PV=""
PV_BATCH_NAME="pvbatch"
AVG_BATCH_NAME="avgbatch"
DEV_BATCH_NAME="devbatch"
PV_BOOL="FALSE"
AVG_BOOL="FALSE"
DEV_BOOL="FALSE"
VARNAME=""
AVGNAME=""
AVGFILE=""
MISSING_FILES="FALSE"

for i in "$@"; do 
  case $i in 
    --ddir=*)
    DATA_DIR="${i#*=}"
    ;;
    --bdir=*)
    BINDIR="${i#*=}"
    ;;
    --h2dir=*)
    H2DIR="${i#*=}"
    ;;
    --h4dir=*)
    H4DIR="${i#*=}"
    ;;
    --list4d=*)
    LIST_4D="${i#*=}"
    ;;
    --list2d=*)
    LIST_2D="${i#*=}"
    ;;
    --listpv=*)
    LIST_PV="${i#*=}"
    ;;
    --pv)
    PV_BOOL="TRUE"
    ;;
    --avg)
    AVG_BOOL="TRUE"
    ;;
    --dev)
    DEV_BOOL="TRUE"
    ;;
    --varname=*)
    VARNAME="${i#*=}"
    ;;
    --avgname=*)
    AVGNAME="${i#*=}"
    ;;
    --avgfile=*)
    AVGFILE="${i#*=}"
    ;;
    --missing)
    MISSING_FILES="TRUE"
    ;;
    *)

    ;;
  esac
done


if [ "$DATA_DIR" == "" ] || [ "$BINDIR" == "" ]; then
  echo "Need to specify data (--ddir) and binary (--bdir) directories"
  exit
fi

if [ ! -e $DATA_DIR ]; then
  mkdir -p $DATA_DIR
fi

if [ "$PV_BOOL" == "FALSE" ] && [ "$AVG_BOOL" == "FALSE" ] && [ "$DEV_BOOL" == "FALSE" ]; then
  echo "Need to specify at least one job flag (--pv, --avg, or --dev)"
  exit
fi

PV_BATCH_NAME=$DATA_DIR/$PV_BATCH_NAME
AVG_BATCH_NAME=$DATA_DIR/$AVG_BATCH_NAME
DEV_BATCH_NAME=$DATA_DIR/$DEV_BATCH_NAME
# This part of the code starts with h2/h4 cam files and generates interpolated variables,
#  then calculates PV that is averaged over 100-500 mb 
if [ "$PV_BOOL" == "TRUE" ]; then 
  
  #check that all directories provided
  if [ "$H4DIR" == "" ] || [ "$H2DIR" == "" ]; then
    echo "Need to provide directories for 2D (--h2dir) and 4D (--h4dir) files" 
    exit
  fi
  #If lists aren't already provided, make them
  if [ "$LIST_4D" == "" ] || [ "$LIST_2D" == "" ] || [ "$LIST_PV" == "" ]; then
    #in each of the directories, make a list of files
    LIST_4D=h4list.txt
    LIST_2D=h2list.txt
    LIST_PV=integ_list.txt

    cd $DATA_DIR
    ls $H2DIR/*_sub.nc > $DATA_DIR/$LIST_2D
    ls $H4DIR/*_sub.nc > $DATA_DIR/$LIST_4D
    #Make a list of the files that will be created
    cd $H4DIR
    ls *_sub.nc | awk -v d="$DATA_DIR" '{sub("_sub.nc",""); print "echo "d"/"$1"_integ.nc"}'| sh > $DATA_DIR/$LIST_PV
  fi

  cd $DATA_DIR
  NUM_H2=$(cat $LIST_2D | wc -l)    
  NUM_H4=$(cat $LIST_4D | wc -l)
   
#  if [ $NUM_H2 -ne $NUM_H4 ]; then
#    echo "Number of h2 and h4 files doesn't match!"
    echo "H4 files: $NUM_H4 H2 files: $NUM_H2"
    echo "Altering list files to only matches."
    MISSING_FILES="TRUE"
    #Check the differences between h2 and h4
    cd $H2DIR
    ls *.nc | rev | cut -d "." -f 2 | rev | cut -d "_" -f 1 > $DATA_DIR/h2_dates.txt
    cd $H4DIR
    ls *.nc | rev | cut -d "." -f 2 | rev | cut -d "_" -f 1 > $DATA_DIR/h4_dates.txt
    cd $DATA_DIR
    comm -1 -2 h2_dates.txt h4_dates.txt > common_dates.txt
    # create new lists with only the common dates
    grep -F -f common_dates.txt $LIST_4D > newh4.txt
    grep -F -f common_dates.txt $LIST_2D > newh2.txt
    grep -F -f common_dates.txt $LIST_PV > newinteg.txt

    LIST_4D=newh4.txt
    LIST_2D=newh2.txt
    LIST_PV=newinteg.txt
#  fi

  cd $DATA_DIR
  #Create a text file with all 3 lists
  paste $LIST_4D $LIST_2D $LIST_PV > combined_list.txt

  #Create the batch script to submit to serial queue
  echo "#!/bin/bash -l" > $PV_BATCH_NAME
  echo "" >> $PV_BATCH_NAME
  echo "#SBATCH -p shared" >> $PV_BATCH_NAME
  echo "#SBATCH -t 20:00:00" >> $PV_BATCH_NAME
  echo "#SBATCH --mem=20GB" >> $PV_BATCH_NAME
  echo "" >> $PV_BATCH_NAME
  echo "cd $DATA_DIR" >> $PV_BATCH_NAME
  echo "cat combined_list.txt | awk '{print \"$BINDIR/BlockingPV --ipl --in \"\$1 \" --in2d \"\$2\" --out \"\$3}'| sh" >> $PV_BATCH_NAME

  if [ "$AVG_BOOL" == "TRUE" ] && [ "$PV_BOOL" == "TRUE" ]; then
    echo "sbatch --dependency=afterok:\${SLURM_JOB_ID} $AVG_BATCH_NAME" >> $PV_BATCH_NAME
  fi
  echo "Submitting PV batch file."
  sbatch $PV_BATCH_NAME
fi

#This part of the code takes the list of input files and generates a single file with the long term daily mean
if [ "$AVG_BOOL" == "TRUE" ]; then
  if [ "$LIST_PV" == "" ]; then 
    echo "Need to specify input list (--listpv=)."
    exit
  fi
  if [ "$AVGNAME" == "" ]; then
    echo "Need to specify average variable name (--avgname=)."
    exit
  fi
  if [ "$VARNAME" == "" ]; then
    echo "Need to specify input variable name (--varname=)."
    exit
  fi
  #The PV list is the input list for the averaging code
  #Create the average file name
  if [ "$AVGFILE" == "" ]; then
    AVG_START=$(head -1 $LIST_PV | cut -d "." -f 1)
    AVGFILE=$AVG_START"_avg_$VARNAME.nc"
    echo "Averaged file name is $AVGFILE"
  fi

  #Create the averaging batch file
  echo "#!/bin/bash -l" > $AVG_BATCH_NAME
  echo "" >> $AVG_BATCH_NAME
  echo "#SBATCH -p shared" >> $AVG_BATCH_NAME
  echo "#SBATCH -t 0:20:00" >> $AVG_BATCH_NAME
  echo "" >> $AVG_BATCH_NAME
  echo "cd $DATA_DIR" >> $AVG_BATCH_NAME
  if [ "$MISSING_FILES" == "TRUE" ]; then
    echo "$BINDIR/BlockingAvg --missing --inlist $LIST_PV --out $AVG_NAME --varname $VARNAME --avgname $AVGNAME" >> $AVG_BATCH_NAME
  else
    echo "$BINDIR/BlockingAvg --missing --inlist $LIST_PV --out $AVGFILE --varname $VARNAME --avgname $AVGNAME" >> $AVG_BATCH_NAME
  fi
  if [ "$DEV_BOOL" == "TRUE" ]; then
    echo "sbatch --dependency=afterok:\${SLURM_JOB_ID} $DEV_BATCH_NAME" >> $AVG_BATCH_NAME
  fi
  if [ "$PV_BOOL" == "FALSE" ]; then
    echo "Submitting averaging batch file."
    sbatch $AVG_BATCH_NAME
  fi
fi

#This part of the code takes the list of input files and the average file and generates deviations from the average
if [ "$DEV_BOOL" == "TRUE" ]; then
  if [ "$LIST_PV" == "" ]; then
    echo "Need to specify input list (--listpv=)."
    exit
  fi
  if [ "$AVGFILE" == "" ]; then
    echo "Need to specify average file (--avgfile=)."
    exit
  fi
  if [ "$VARNAME" == "" ]; then
    echo "Need to specify input variable name (--varname=)."
    exit
  fi
  if [ "$AVGNAME" == "" ]; then
    echo "Need to specify average variable name (--avgname=)."
    exit
  fi
  #Create the deviations batch file
  echo "#!/bin/bash -l" > $DEV_BATCH_NAME
  echo "" >> $DEV_BATCH_NAME
  echo "#SBATCH -p shared" >> $DEV_BATCH_NAME
  echo "#SBATCH -t 2:00:00" >> $DEV_BATCH_NAME
  echo "" >> $DEV_BATCH_NAME
  echo "cd $DATA_DIR" >> $DEV_BATCH_NAME
  echo "$BINDIR/BlockingDevs --inlist $LIST_PV --avg $AVGFILE --varname $VARNAME --avgname $AVGNAME" >> $DEV_BATCH_NAME

  if [ "$AVG_BOOL" == "FALSE" ]; then 
    echo "Submitting deviations batch file."
    sbatch $DEV_BATCH_NAME
  fi

fi
