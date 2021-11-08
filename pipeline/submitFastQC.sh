#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input dir> <output dir>"
    exit
fi

if [ ! -d $1 ]; then
    echo "error: input directory not found: $1" >&2
    exit
fi

if [ ! -d $2 ]; then
    echo "error: output directory not found: $2" >&2
    exit
fi

IN=$1
OUT=$2

for f in `find "$IN" -name "*.f*q.gz"`; do
    bash ../../../pipeline/runFastQC.sh "$OUT" "$f"
done
