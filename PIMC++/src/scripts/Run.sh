#!/bin/sh
echo "Running" $1 "with" $3 "from" $2 "..."
$1 $3 1>>$3.out 2>&1
