#!/bin/bash

# GH_script.sh
# Version 1.0 Jan 19, 2017
# Author: Marielle Pinheiro

# This script automates the interpolation/extraction
# of the Z500 variable from a model level geopotential height
# variable (cam5 runs)

# TO RUN THE SCRIPT:
# At the command line enter:
# ./GH_script.sh WITH THE ADDITIONAL FLAGS:
# --ddir=[PATH WHERE DATA WILL BE STORED] ##MANDATORY
# --bdir=[PATH WHERE BINARIES ARE LOCATED] ##MANDATORY
# JOB FLAGS:
# --Z500: run the binary Var4Dto3D (extracts 500 mb height)
# --avg: find long term daily average using provided list of Z500 files
# --dev: find the deviations from the long term daily average
# FOR Z500 CALCULATIONS ONLY:
# --h2dir=[PATH WHERE H2 FILES LOCATED] ##MANDATORY
# --h4dir=[PATH WHERE H4 FILES LOCATED] ##MANDATORY
# FOR AVERAGE CALCULATIONS:
# --listz=[LIST OF INPUT FILES] ##MANDATORY IF NOT DOING Z500
# --varname=[NAME OF VARIABLE BEING AVERAGED] ##MANDATORY
# --avgfile=[NAME OF AVERAGING FILE] ##MANDATORY
# --missing: specify if there are gaps in the data (otherwise
#         the averaging code will prematurely terminate!)
# FOR THE DEVIATION CALCULATIONS:
# --listz=[LIST OF INPUT FILES] ##MANDATORY IF NOT DOING Z500
# --avgfile=[NAME OF AVERAGED FILE] ##MANDATORY IF NOT RUNNING AVG

DATA_DIR=""
BINDIR=""
H2DIR=""
H4DIR=""
LIST_4D=""
LIST_2D=""
LIST_Z=""
LIST_BLOBS=""
Z_BATCH_NAME="zbatch"
AVG_BATCH_NAME="avgzbatch"
DEV_BATCH_NAME="devzbatch"
BLOBS_BATCH_NAME="blobszbatch"
Z_BOOL="FALSE"
AVG_BOOL="FALSE"
DEV_BOOL="FALSE"
BLOBS_BOOL="FALSE"
VARNAME=""
AVGNAME=""
DEVNAME=""
BLOBSNAME=""
AVGFILE=""
BLOBSFILE=""

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
    --listz=*)
    LIST_Z="${i#*=}"
    ;;
    --listblobs=*)
    LIST_BLOBS="${i#*=}"
    ;;
    --Z500)
    Z_BOOL="TRUE"
    ;;
    --avg)
    AVG_BOOL="TRUE"
    ;;
    --dev)
    DEV_BOOL="TRUE"
    ;;
    --blobs)
    BLOBS_BOOL="TRUE"
    ;;
    --varname=*)
    VARNAME="${i#*=}"
    ;;
    --avgname=*)
    AVGNAME="${i#*=}"
    ;;
    --devname=*)
    DEVNAME="${i#*=}"
    ;;
    --blobsname=*)
    BLOBSNAME="${i#*=}"
    ;;
    --avgfile=*)
    AVGFILE="${i#*=}"
    ;;
    --blobsfile=*)
    BLOBSFILE="${i#*=}"
    ;;
    *)

    ;;
  esac
done

if [ "$VARNAME" == "" ]; then
  VARNAME="Z3"
fi
if [ "$AVGNAME" == "" ]; then
  AVGNAME="AVGZ"
fi


if [ "$DATA_DIR" == "" ] || [ "$BINDIR" == "" ]; then
  echo "Need to specify data (--ddir) and binary (--bdir) directories"
  exit
fi

if [ ! -e $DATA_DIR ]; then
  mkdir -p $DATA_DIR
fi

if [ "$Z_BOOL" == "FALSE" ] && [ "$AVG_BOOL" == "FALSE" ] && [ "$DEV_BOOL" == "FALSE" ]; then
  echo "Need to specify at least one job flag (--Z500, --avg, or --dev)"
  exit
fi

Z_BATCH_NAME=$DATA_DIR/$Z_BATCH_NAME
AVG_BATCH_NAME=$DATA_DIR/$AVG_BATCH_NAME
DEV_BATCH_NAME=$DATA_DIR/$DEV_BATCH_NAME


#This part of the code creates the Z500 files
if [ "$Z_BOOL" == "TRUE" ]; then
  #Check that all the directories have been provided
  if [ "$H4DIR" == "" ] || [ "$H2DIR" == "" ]; then
    echo "Need to provide h2 and h4 directory locations."
    exit
  fi
  #If the lists aren't already provided, then make them
  if [ "$LIST_4D" == "" ] || [ "$LIST_2D" == "" ] || [ "$LIST_Z" == "" ]; then
    LIST_4D=h4_zlist.txt
    LIST_2D=h2_zlist.txt
    LIST_Z=z500_list.txt

    cd $H2DIR
    ls *_sub.nc | rev | cut -d "." -f 2 | rev | cut -d "_" -f 1 > $DATA_DIR/h2_zdates.txt
    cd $H4DIR
    ls *_sub.nc | rev | cut -d "." -f 2 | rev | cut -d "_" -f 1 > $DATA_DIR/h4_zdates.txt
    ls *_sub.nc | awk -v d="$DATA_DIR" '{sub("_sub.nc",""); print "echo "d"/"$1"_z500.nc"}' | sh > $DATA_DIR/$LIST_Z
    cd $DATA_DIR
    comm -1 -2 h2_zdates.txt h2_zdates.txt > common_zdates.txt

    ls $H2DIR/*_sub.nc > $DATA_DIR/$LIST_2D
    ls $H4DIR/*_sub.nc > $DATA_DIR/$LIST_4D
    
    #Make a list with only the common dates
    grep -F -f common_zdates.txt $LIST_4D > new_zh4.txt
    grep -F -f common_zdates.txt $LIST_2D > new_zh2.txt
    grep -F -f common_zdates.txt $LIST_Z > new_z500.txt
  fi

  #Create a text file containing all 3 lists
  paste $LIST_4D $LIST_2D $LIST_Z > combined_zlist.txt
  echo "#!/bin/bash -l" > $Z_BATCH_NAME
  echo "" >> $Z_BATCH_NAME
  echo "#SBATCH -p shared" >> $Z_BATCH_NAME
  echo "#SBATCH -o zbatch.output" >> $Z_BATCH_NAME
  echo "#SBATCH -t 2:00:00" >> $Z_BATCH_NAME
  echo "#SBATCH -C haswell" >> $Z_BATCH_NAME
  echo "#SBATCH -L SCRATCH" >> $Z_BATCH_NAME
  echo "" >> $Z_BATCH_NAME
  echo "cd $DATA_DIR" >> $Z_BATCH_NAME
  echo "cat combined_zlist.txt | awk '{print \"$BINDIR/Var4Dto3D --ipl --varlist $VARNAME --in \"\$1 \" --in2d \"\$2\" --out \"\$3}'| sh" >> $Z_BATCH_NAME

  if [ "$AVG_BOOL" == "TRUE" ] && [ "$Z_BOOL" == "TRUE" ]; then
    echo "sbatch --dependency=afterok:\${SLURM_JOB_ID} $AVG_BATCH_NAME" >> $Z_BATCH_NAME
  fi
  echo "Submitting Z500 batch file."
#  sbatch $Z_BATCH_NAME

fi

#This part of the code takes the list of input files and generates a single file with the long term daily mean

if [ "$AVG_BOOL" == "TRUE" ]; then
  if [ "$LIST_Z" == "" ]; then
    echo "Need to specify input list."
    exit
  fi
  #The Z list is the input list for the averaging code
  #Create the averaged file name
  if [ "$AVGFILE" == "" ]; then
    AVG_START=$(head -1 $LIST_Z | cut -d "." -f 1,2)
    AVG_FILE=$AVG_START".Z500_avg.nc"
  fi
  cd $DATA_DIR
  #Create the batch file
  echo "#!/bin/bash -l" > $AVG_BATCH_NAME
  echo "" >> $AVG_BATCH_NAME
  echo "#SBATCH -p shared" >> $AVG_BATCH_NAME
  echo "#SBATCH -o avg_batch.output">> $AVG_BATCH_NAME
  echo "#SBATCH -t 0:30:00" >> $AVG_BATCH_NAME
  echo "#SBATCH -C haswell" >> $AVG_BATCH_NAME
  echo "#SBATCH -L SCRATCH" >> $AVG_BATCH_NAME
  echo "" >> $AVG_BATCH_NAME
  echo "cd $DATA_DIR" >> $AVG_BATCH_NAME
  echo "$BINDIR/BlockingAvg --missing --inlist $LIST_Z --out $AVGFILE --varname $VARNAME --avgname $AVGNAME" >> $AVG_BATCH_NAME

  if [ "$DEV_BOOL" == "TRUE" ]; then
    echo "sbatch --dependency=afterok:\${SLURM_JOB_ID} $DEV_BATCH_NAME" >> $AVG_BATCH_NAME
  fi
  if [ "$Z_BOOL" == "FALSE" ]; then
  echo "Submitting averaging batch file."
#  sbatch $AVG_BATCH_NAME
  fi
fi

#This part of the code takes the list of input files and the average file and generates deviations from the average
if [ "$DEV_BOOL" == "TRUE" ]; then
  if [ "$LIST_Z" == "" ]; then
    echo "Need to specify input list."
    exit
  fi
  if [ "$AVGFILE" == "" ]; then
    echo "Need to specify average file."
    exit
  fi
  cd $DATA_DIR
  #Create the deviations batch file
  echo "#!/bin/bash -l" > $DEV_BATCH_NAME
  echo "" >> $DEV_BATCH_NAME
  echo "#SBATCH -p shared" >> $DEV_BATCH_NAME
  echo "#SBATCH -o dev_batch.output">>$DEV_BATCH_NAME
  echo "#SBATCH -t 2:00:00" >> $DEV_BATCH_NAME
  echo "#SBATCH -C haswell">> $DEV_BATCH_NAME
  echo "#SBATCH -L SCRATCH" >> $DEV_BATCH_NAME
  echo "" >> $DEV_BATCH_NAME
  echo "cd $DATA_DIR" >> $DEV_BATCH_NAME
  echo "$BINDIR/BlockingDevs --gh --inlist $LIST_Z --avg $AVGFILE --varname $VARNAME --avgname $AVGNAME" >> $DEV_BATCH_NAME
  if [ "$AVG_BOOL" == "FALSE" ]; then
    echo "Submitting deviations batch file."
#    sbatch $DEV_BATCH_NAME
  fi

fi



