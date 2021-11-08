#!/bin/bash

set -e

IN="/mnt/picea/storage/projects/08_ABA/data/sortmerna"
OUT="/mnt/picea/storage/projects/08_ABA/data/trimmomatic"

trim="SLIDINGWINDOW:5:30 MINLEN:30"

export SLURM_SUBMIT_DIR="."

for f in `find "$IN" -name "*.f*q.gz"`; do
    fname=${f##*/}
    time bash ../../../pipeline/runTrimmomatic.sh -s "$f" "$OUT" 33 $trim > "$OUT/$fname.log" 2>&1
done