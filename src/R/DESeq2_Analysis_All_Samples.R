#'---
#' title: "DESeq2 Analysis of All Samples"
#' author: "Crystal Flores"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---


library(DESeq2)
#'Location of the HTSeq files provided by Nico 
directory  <- "/mnt/picea/projects/spruce/29_Spruce_Seeds_Project/htseq-nico"

dir(directory)

#' these 2 commands do the same 
#sampleFiles <- list.files(directory)[grep("STAR", list.files(directory))]
sampleFiles <- grep("STAR", list.files(directory),value = TRUE)
sampleName <- sub(".*_P1790_","",sub("_sortmerna.*\\.txt","",list.files(directory)))
sampleCsv <- read.csv("~/Git/UPSCb/projects/spruce-seeds/doc/samples.csv")
sampleCsv$sampleIndex <- as.character(seq(101,112))
sampleCondition <- sampleCsv[match(sampleName,sampleCsv$sampleIndex),"Tissue"]
#'To check what type of value we have we can command class()
class(sampleCondition)

#'If the class were a factor it looks like this
#as.integer(sampleCondition)
#'If the class were a character it would look like this 
#as.character(sampleCondition)

sampleTable <- data.frame(sampleName=sampleName,
                          fileName=sampleFiles,
                          condition=sampleCondition)                       
sampleTable

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                       directory= directory,
                                       design= ~ condition)
ddsHTSeq


#'Differential Expression Analysis 

dds<-DESeq(ddsHTSeq)
res <- results(dds)
res

#'Order results by smallest adjusted p-value
resOrdered<- res[order(res$padj),]
head(resOrdered)

#'Result Summary 
summary(res)

#'MA Plot 
resMLE<-results(dds,addMLE=TRUE)
head(resMLE,4)
plotMA(resMLE,main="DESeq2", ylim=c(-2,2))
#'Making the plot interactive in order to detect the row no. of individual genes
#'remember row no. is the genes that we are looking at, columns are our 12 samples
#identify(res$baseMean, res$log2FoldChange)

#'Plot Counts
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition",returnData=TRUE)
d
library("ggplot2")
ggplot(d,aes(x=condition, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(100,200,20000,30000))

#'Information represented in the columns of the results
mcols(res)$description

head(results(dds, addMLE=TRUE),4)

