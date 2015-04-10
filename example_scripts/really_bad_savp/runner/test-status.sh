#!/bin/bash

FILENAME=$1

if [ -f "status-$FILENAME" ];
then
    rm -rf "status-$FILENAME"
fi

DATE1=$(date +%s)
while read p; do
    NAME=$(echo $p | cut -f 1 -d " " )
    JOB=$(echo $p | cut -f 2 -d " " )
    echo "name: $NAME"
    STATUS=$(jobs-status $JOB | cut -f 1)
    echo -e "$NAME\t$JOB\t$STATUS" >> "status-$FILENAME"
done <$FILENAME

DATE2=$(date +%s)

ESEC=`expr $DATE2 - $DATE1`
echo "$ESEC seconds elapsed running this script"
