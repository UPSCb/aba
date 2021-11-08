#' ---
#' title: "Extract P.tremula_X_Ptremuloides synthetic transcripts"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' ## Environment
#' Set the working dir
setwd("/mnt/picea/storage/reference/Populus-tremula_X_Populus-tremuloides/v1.0/fasta")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/storage/reference/Populus-tremula_X_Populus-tremuloides/v1.0/fasta")
#' ```
#' Load libraries
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(genomeIntervals))

#' Source helper
source("~/Git/UPSCb/src/R/sequenceUtilities.R")

#' # Extract
#' Read the genome
genome <- readDNAStringSet("Potrx01-genome.fa")

#' Read the gff3 file
gff3 <- readGff3("../gff3/Potrx01-synthetic-transcripts.gff3.gz")

#' Get all transcripts
trx <- extractMrnaFromGenome(gff3,genome)

#' # Export
writeXStringSet(trx,file="Potrx01-synthetic-transcripts.fa")

#' # Session Info
sessionInfo()
