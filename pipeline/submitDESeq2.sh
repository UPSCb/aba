#!/bin/bash

set -e

htseqdir="/mnt/picea/storage/projects/08_ABA/htseq"
IN="/mnt/picea/storage/projects/08_ABA/analysis/deseq/metadata"
OUT="/mnt/picea/storage/projects/08_ABA/analysis/deseq2"

for f in `find "$IN" -name "*.txt"`; do
    echo "Running DESeq2 on ${f##*/}"
    Rscript ../src/runDESeq2.R "$f" "$htseqdir" "$OUT"
done
