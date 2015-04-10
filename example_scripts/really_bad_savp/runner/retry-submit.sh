#!/bin/bash

JOBLISTFILE=$1
while read p; do
    JOB=$(echo $p | sed -e 's/^[[:space:]]*//g' -e 's/[[:space:]]*\$//g' )
    echo $JOB
done <$JOBLISTFILE

