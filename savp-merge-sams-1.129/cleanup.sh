#!/bin/bash

# This is a BAM file you wish to exclude from deletion
outputFileName=$1

# Cleverly find and delete any BAM or SAM that is NOT the output file
for FF in $(find -regextype awk -iregex ".*/*.bam|.*/*.sam" | grep -v $outputFileName)
do
    RMF=$(basename $FF)
    echo "$Detected $RMF as input file."
    rm -rf $RMF
done
