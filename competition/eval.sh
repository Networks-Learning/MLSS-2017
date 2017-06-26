#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

NETWORK_FILE=$1
SOLUTION_FILE=$2
PROB=0.5

sources=$(cat "$SOLUTION_FILE")

op=$(../code/mlss-17/simulate_ic -c:5000 -sg:0 -st:0 -sm:0 -v:0 -p:"$PROB" -t:0 -s:"$sources" -o:"/tmp/cascades" -i:"$NETWORK_FILE" 2>/dev/null | fgrep 'Average')

echo $op
