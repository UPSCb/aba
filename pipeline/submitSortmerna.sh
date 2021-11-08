#!/bin/bash

if [ $# != 1 ]; then
    echo "Usage: $0 input_dir" > /dev/stderr
    exit 1
fi

if [ ! -d $1 ]; then
    echo "First argument must be an existing directory" > /dev/stderr
    exit 2
fi

# Make sure the merge-paired-reads.sh is in the path
if [[ ! `which merge-paired-reads.sh` ]]; then
    export PATH="/mnt/picea/home/niklasm/git/UPSCb/src/sh:$PATH"
fi

export SORTMERNADIR="../../../data/sortmerna"

IN=$1

TMP="/mnt/picea/storage/projects/08_ABA/tmp/sortmerna"
OUT="/mnt/picea/storage/projects/08_ABA/data/sortmerna"

if [ ! -e ../../../data/sortmerna ]; then
    echo "ERROR sortmerna data not found, link (or copy) it!" > /dev/stderr
    exit 1
fi

for sample in `find "$IN" -name "*.f*q.gz"`; do
    echo "Running sortmerna for $sample..."
    bash ../../../pipeline/runSortmerna.sh -u "$OUT" "$TMP" "$sample"
done
