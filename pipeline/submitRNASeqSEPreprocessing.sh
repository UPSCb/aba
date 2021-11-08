#!/bin/bash -l

set -ex

proj=u2015021

mail="katja.stojkovic@slu.se"

in=/mnt/picea/projects/aspseq/rbhalerao/ABA/Potra01/raw

out=/mnt/picea/projects/aspseq/rbhalerao/ABA/Potra01

genome=/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/STAR/2.5.2b/Potra01

gff3=/mnt/picea/storage/reference/Populus-tremula/v1.1/gff3/Potra01-gene-synthetic-transcripts-wo-intron.gff3

gtf=/mnt/picea/storage/reference/Populus-tremula/v1.1/gtf/Potra01-gene-mRNA-wo-intron.gtf

kallistoFasta=/mnt/picea/storage/reference/Populus-tremula/v1.1/fasta/Potra01-mRNA.fa
			  
kallistoIndex=/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/kallisto/Potra01-mRNA.fa.inx

start=5
end=9
mem=256

module load bioinfo-tools FastQC Trimmomatic sortmerna star samtools htseq fastQvalidator kallisto

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

# if you need strand specificity add -t before $proj
for f in `find $in -name "*.fastq.gz"`; 
do
  bash $UPSCb/pipeline/runRNASeqSEPreprocessing.sh -s $start -e $end -m $mem -g $genome \
  -G $gtf -H $gff3 -f $kallistoFasta -K $kallistoIndex -T "SLIDINGWINDOW:5:20 MINLEN:35ob" -M 200 -S 20 $proj $mail $f $out
done