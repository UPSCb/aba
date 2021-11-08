#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Input folder required"
    exit
fi

IN=$1

for f in `find "$IN" -name "*.f*q.gz"`; do
    fastQValidator --file "$f" --maxErrors 100 --printableErrors 100
    if [ $? -ne 0 ]; then
        echo "$f is not a valid FastQ file"
    fi
done
