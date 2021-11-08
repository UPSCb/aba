#!/bin/bash

set -e

htseqdir="/mnt/picea/storage/projects/08_ABA/htseq"
IN="/mnt/picea/storage/projects/08_ABA/analysis/deseq/metadata"
OUT="/mnt/picea/storage/projects/08_ABA/analysis/deseq"

for f in `find "$IN" -name "*.txt"`; do
    echo "Running DESeq on ${f##*/}"
    Rscript ../src/runDESeq.R "$f" "$htseqdir" "$OUT"
done
