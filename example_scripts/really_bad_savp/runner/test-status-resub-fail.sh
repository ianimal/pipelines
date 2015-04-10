#!/bin/bash

# File derived from test-submit
FILENAME=$1
# Directory containing JOBNAME.json files
JOBDIR=$2
# Stem for new submit list
JOBLISTNAME=$3

DATE1=$(date +%s)
while read p; do
    NAME=$(echo $p | cut -f 1 -d " " )
    JOB=$(echo $p | cut -f 2 -d " " )
    echo "name: $NAME"
    STATUS=$(jobs-status $JOB | cut -f 1)
    if [[ $STATUS =~ "FAILED" ]];
    then
    # Resubmit
        echo "Submitting $NAME again"
        JOB=""
        RETRY=1
        while [[ ! $JOB =~ "0001-007" ]]
        do
            echo "...Attempt $RETRY"
            JOB=$( jobs-submit -F "${JOBDIR}/${NAME}.json" | cut -f 4 -d " " -)
            let RETRY=RETRY+1
        done
        echo -e "$NAME\t$JOB" >> "submit-jobs-${JOBLISTNAME}.txt"

    fi

done <$FILENAME

DATE2=$(date +%s)

ESEC=`expr $DATE2 - $DATE1`
echo "$ESEC seconds elapsed running this script"
