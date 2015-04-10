#!/bin/bash

JOBDIR=$1
JOBLISTNAME=$2

rm -rf ${JOBLISTNAME}
COUNTER=0
DATE1=$(date +%s)
for FF in $(ls $JOBDIR)
do
    echo "Submitting job $FF"
    JNAME="${FF%%.*}"

    JOB=""
    RETRY=1
    while [[ ! $JOB =~ "0001-007" ]]
    do
        echo "...Attempt $RETRY"
        JOB=$( jobs-submit -F $JOBDIR/$FF | cut -f 4 -d " " -)
        let RETRY=RETRY+1
    done

    echo -e "$JNAME\t$JOB" >> "submit-jobs-${JOBLISTNAME}.txt"
    let COUNTER=COUNTER+1
done
DATE2=$(date +%s)

ESEC=`expr $DATE2 - $DATE1`
echo "$ESEC seconds elapsed submitting $COUNTER jobs"

exit 0


