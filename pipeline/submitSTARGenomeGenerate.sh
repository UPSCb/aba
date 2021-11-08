#!/bin/bash -l

## error and verbose
set -ex

## define a function
usage () {
echo >&2 \
"This function take one mandatory argument as parameter:
     the index to create; one of 'Potra01', 'Potrs01', 'Potri03' or 'Potrx01'
"
exit 1
}

## args number
if [ $# != 1 ]; then
    usage
    exit 1
fi

case "$1" in
    Potra01)
	inxDir=/mnt/picea/storage/reference/Populus-tremula/v1.0/indices/STAR
	fasta=/mnt/picea/storage/reference/Populus-tremula/v1.0/fasta/Potra01-genome.fa
	gtf=/mnt/picea/storage/reference/Populus-tremula/v1.0/gtf/Potra01-gene-wo-intron.gtf
	;;
    Potrs01)
	inxDir=/mnt/picea/storage/reference/Populus-tremuloides/v1.0/indices/STAR
	fasta=/mnt/picea/storage/reference/Populus-tremuloides/v1.0/fasta/Potrs01-genome.fa
	gtf=/mnt/picea/storage/reference/Populus-tremuloides/v1.0/gtf/Potrs01-gene-wo-intron.gtf
	;;
    Potri03)
	inxDir=/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/indices/STAR
	fasta=/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/fasta/Ptrichocarpa_v3.0_210.fa
	gtf=/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/gtf/Ptrichocarpa_v3.0_210_gene_exons.gtf
	;;
	  Potrx01)
	inxDir=/mnt/picea/storage/reference/Populus-tremula_X_Populus-tremuloides/v1.0/indices/STAR
	fasta=/mnt/picea/storage/reference/Populus-tremula_X_Populus-tremuloides/v1.0/fasta/Potrx01-genome.fa
	gtf=/mnt/picea/storage/reference/Populus-tremula_X_Populus-tremuloides/v1.0/gtf/Potrx01-genome.gtf
	;;
    *)
	usage;;
esac

## prepare
mkdir -p $inxDir/$1
sbatch -e $inxDir/$1.err -o $inxDir/$1.out --mem=200G -n 32 --mail-user="nicolas.delhomme@umu.se" $UPSCb/pipeline/runSTARGenomeGenerate.sh -m 200000000000 -p 32 -f gtf $inxDir/$1 $fasta $gtf
