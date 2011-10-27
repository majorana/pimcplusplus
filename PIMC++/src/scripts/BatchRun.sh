#!/bin/sh

echo "Using" $1 "Processes"
nLines=$(wc -l < $3 )
nInputLines=$(($nLines-1))
echo $3 "has" $nInputLines "input lines"

cat $3 | xargs -n1 -P$1 sh ${0%/*}/Run.sh $2 $3

# mutt -s "Batch Run Completed" ethan.w.brown@gmail.com < $3
