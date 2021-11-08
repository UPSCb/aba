#!/bin/bash

## stop on error and be verbose
set -ex

## global vars
proj=u2015020
mail="katja.stojkovic@slu.se"

# a usage function
usage () {
echo >&2 \
"This function take one mandatory argument as parameter:
     the index to create; one of 'Potrx01', 'Potri03'
"
exit 1
}

# check the arguments
if [ $# != 1 ]; then
    usage
fi

# use the argument
case "$1" in
  Potri03)
    GFF="/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/gff/Ptrichocarpa_v3.0_210_synthetic-gene-models.gff3"
    in=/mnt/picea/projects/aspseq/ABA/star
    out=/mnt/picea/projects/aspseq/ABA/htseq
  ;;
  Potrx01)
    in=/mnt/picea/projects/aspseq/rbhalerao/ABA/Potrx01/star
    out=/mnt/picea/projects/aspseq/rbhalerao/ABA/Potrx01/analysis/HTSeq_T89andABI1/
    GFF=/mnt/picea/projects/aspseq/rbhalerao/histone-methylation/Reference/Potrx01-gene_16-11-18_corrected.gff3
  ;;
  ?)
  echo "unknown option: $1"
  usage
  ;;
esac

# check the out dir
if [ ! -d $out ]; then
  mkdir -p $out
fi

# check the UPSCb env var
if [ -z $UPSCb ]; then
  echo "The UPSCb env. var. is not defined"
  usage
fi

# process
for f in `find $in -name "*_sortmerna_trimmomatic_STAR.cram"`; do
  fnam=`basename ${f//.cram/}`
    sbatch --mail-type=ALL --mail-user $mail -o $out/$fnam.out -e$out/$fnam.err -J S-$fnam \
    $sbatchOpt $UPSCb/pipeline/runHTSeq.sh $out $f $GFF
    
done

