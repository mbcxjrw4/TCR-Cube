#!/bin/bash

# runAntigenIdentification.sh
# single machine version of APP
# author: Ruoxi Wang - ruoxi.wang@adaptimmune.com
#
# a txt mandate file that looks like ...
# Indication
# adrenocortical_cancer
# 
# usage
# review params below, then
# ./runAntigenIdentification 3
# to run the script usng the parameters listed on the third line of the mandate file
# or
# seq 1 180 | xargs -P 2 -n l bash runAntigenIdentification.sh
# where n is the number of lines in the MANDATE file
#       P is the maximum number of processors to use at once
# below, specify a full path to the input data directory
# and a mandate file that lists all of the sample names
# and a full path to a directory where the results will be placed
#

MANDATE="/data/input/indication_target.txt"
OUTPUT_DIR="data/processed/intermediate/"

###
# no changes required beyond this point
###

APP_NAME=$0
TASK_ID=$1

#retrieve the variables for the job from the mandate file
THIS_MANDATE=$(sed -n -e "${TASK_ID} p" ${MANDATE})

# directory named by indication for intermediate files
mkdir -p ${THIS_MANDATE}

#copy the Rmd file to temporary directory and enter into it
cp identify_pairs.R ${THIS_MANDATE}/
cd ${THIS_MANDATE}

#log what will be done for this specfic task
echo -e "starting ${APP_NAME} for ${THIS_MANDATE}." \

START=$(date +%s)

# call R script for $SUBJECT_ID $STUDY_ID
EXIT_CODE=$(Rscript identify_pairs.R ${THIS_MANDATE}) 

# remove the Rmd file and leave the directory
rm -r *.R
cd ..

# move the result to data/processed/intermediate/
mv -rf ${THIS_MANDATE}/ ${OUTPUT_DIR}/

#log what happened
END=$(date +%s)
TOTAL_TIME=$(( $END - $START )) # $(echo "$END - $START" | bc)
echo "exited with code ${EXIT_CODE} for ${THIS_MANDATE} after ${TOTAL_TIME} seconds total time" 

exit 0
