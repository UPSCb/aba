#!/bin/sh
set -ex

## create the dirs
cd /mnt/picea/projects/aspseq/ABA
mkdir Potrx01
ln -s ../data/trimmomatic
mkdir star
mkdir htseq
