#!/bin/bash

FILENAME=$1

DATE1=$(date +%s)

while read p; do
    JOB=$(echo $p | sed -e 's/^[[:space:]]*//g' -e 's/[[:space:]]*\$//g' )
    echo "job: $JOB"
    jobs-delete $JOB
done <$FILENAME

DATE2=$(date +%s)

ESEC=`expr $DATE2 - $DATE1`
echo "$ESEC seconds elapsed running this script"
