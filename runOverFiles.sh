#!/bin/bash

if [ "$2" = "" ]
then 
    echo "usage `basename $0` <dir name> <file base>"
    exit 1
fi

DIR_NAME=$1
BASE_NAME=$2



#for startVal in 5001 6001 7001 8001 9001; do
for startVal in `seq 1 1000 10000`; do
    let endVal=$startVal+999
    echo $startVal $endVal
    OUTFILE=${DIR_NAME}/pca/pca_${BASE_NAME}_${startVal}_${endVal}.root
    time ./makePCAFile $DIR_NAME $BASE_NAME 1000 $OUTFILE $startVal > log.txt 2>&1 
done


#./makePCAFile /unix/anita1/creamtea/container13m_10cmtarget/ container_target 1000 /unix/anita1/creamtea/container13m_10cmtarget/pca/pca_container_8001_9000files.root 8001 > log.txt
