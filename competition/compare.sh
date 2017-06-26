#!/bin/bash

SRC=$1
REF=$2

files=("network-130825-1"
       "network-130825-10"
       "network-130825-2"
       "network-130825-3"
       "network-130825-4"
       "network-130825-5"
       "network-130825-6"
       "network-130825-7"
       "network-130825-8"
       "network-130825-9"
       "network-130826-1"
       "network-130826-10"
       "network-130826-2"
       "network-130826-3"
       "network-130826-4"
       "network-130826-5"
       "network-130826-6"
       "network-130826-7"
       "network-130826-8"
       "network-130826-9"
       "network-130827-1"
       "network-130827-10"
       "network-130827-2"
       "network-130827-3"
       "network-130827-4"
       "network-130827-5"
       "network-130827-6"
       "network-130827-7"
       "network-130827-8"
       "network-130827-9"
       "network-130828-1"
       "network-130828-10"
       "network-130828-2"
       "network-130828-3"
       "network-130828-4"
       "network-130828-5"
       "network-130828-6"
       "network-130828-7"
       "network-130828-8"
       "network-130828-9")

for ii in "${files[@]}"; do
  for s in $SRC/$ii-*; do
    for d in $REF/$ii-*; do
      echo "$s" "$d";
      NETWORK_FILE=$(echo $ii | sed -e 's/-[0-9]*$//')
      ./eval.sh "../code/mlss-17/data/${NETWORK_FILE}.txt" $s
      ./eval.sh "../code/mlss-17/data/${NETWORK_FILE}.txt" $d
      # diff -u "$s" "$d";
      echo "************************";
      echo;
    done;
  done;
done;

