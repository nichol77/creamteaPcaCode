#!/bin/bash

if [ "$2" = "" ]
then 
    echo "usage `basename $0` <dir name> <file base>"
    exit 1
fi

DIR_NAME=$1
BASE_NAME=$2



#for startVal in 5001 6001 7001 8001 9001; do
for indVal in `seq 1 10`; do
    let startVal=$indVal-1
    let startVal=100*$startVal
    let startVal=$startVal+1
    let endVal=$startVal+99
#    echo $startVal $endVal
    OUTFILE=${DIR_NAME}/pca/pca_${BASE_NAME}_million_$indVal.root
    echo $OUTFILE
    time ./makePCAFile $DIR_NAME $BASE_NAME 100 $OUTFILE $startVal > log_${BASE_NAME}.txt 2>&1 

done


#./makePCAFile /unix/anita1/creamtea/container13m_10cmtarget/ container_target 1000 /unix/anita1/creamtea/container13m_10cmtarget/pca/pca_container_8001_9000files.root 8001 > log.txt
