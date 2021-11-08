#' ---
#' title: "Quality assessment of the RNA-seq samples from ABA experiment"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/aspseq/rbhalerao/ABA/Potrx01/analysis/all_samples/HTSeq")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/aspseq/rbhalerao/ABA/Potrx01/analysis/HTSeq")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(stringr))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Create a palette
pal <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(12,"Set3")

#' Register the default plot margin
mar <- par("mar")

#' Construct sample information
files <- dir(pattern = ".cram.txt")
samples <- data.frame(SampleName = str_extract(files, "^[0-9]+w_[a-z]+[0-9]+_[0-9]_L[0-9]"),
                      Seq_laneID = str_extract(files, "L[0-9]"),
                      plantline = str_extract(files, "abi1|t89"),
                      time_point = str_extract(files, "^[0-9]+w"),
                      BioRep = str_extract(files, "^[0-9]+w_[a-z]+[0-9]+_[0-9]"))
# write the sample info to a file
write.csv(samples, file = "~/Git/UPSCb/projects/aba/doc/sampleInfo_allSamples_withTechRep.csv")
  
#' Read the HTSeq files in a matrix separatly from the two experiment 
res <- mclapply(dir(pattern="\\.txt",full.names=TRUE),function(f){
  read.delim(f,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=max(mcaffinity()))

names(res) <- sub("_sortmerna.*\\.txt","",dir(pattern="\\.txt"))

#' Reorder the sample data.frame according to the way
#' the results were read in
samples <- samples[match(names(res),samples$SampleName),]

#' Raw Data QC analysis
addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]

#' Extract the HTSeq stat lines
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

#' Convert them into percentages
apply(count.stats,2,function(co){round(co*100/sum(co))})

#' Plot the stats
col <- pal[1:nrow(count.stats)]
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=col,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,4e+6),cex.names=.6)
legend("top",fill=col,legend=gsub("_"," ",rownames(count.stats)),bty="n",cex=0.8)
par(mar=mar)

#' The average percentage of aligned reads is %
boxplot(unlist(count.stats["aligned",]/colSums(count.stats)),
        main="aligned reads",ylab="percent aligned")

#' Check how many genes are never expressed
sel <- rowSums(count.table) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(count.table),digits=1),
        sum(sel),
        nrow(count.table))

#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative coverage is as expected
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples coloured by time point

plot.multidensity(log10(count.table),
                  col=pal[as.integer(factor(samples$time_point))],
                  legend.x="topright",
                  legend.col=pal[as.integer(factor(samples$time_point))],
                  legend.lwd=1,
                 legend.cex = 0.4,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' Remove genes that are never expressed and save the raw counts
count.table <- count.table[rowSums(count.table)>0, ]
write.csv(count.table,"raw-unnormalised-data-t89aba_withTechRep.csv")

#' Join technical replicates  
#' a) in count table
count.table <- do.call(cbind,
  lapply(split.data.frame(t(count.table), samples$BioRep), colSums))

write.csv(count.table, file = "raw-unnormalised-data-t89aba_bioRep.csv")

#' b) in sample information
samples <- samples[, -c(1,2)]
samples <- unique(samples)
write.csv(samples, file= "~/Git/UPSCb/projects/aba/doc/sampleInfo_allSamples_BioRep.csv")


#' ## QC 

#' # Data normalisation 
#'  For visualization, the data is
#' submitted to a variance stabilization
#' transformation using DESeq2. The 
#' dispersion is estimated independently
#' of the sample time points and replicate

#' Create groups representing condition
samples$condition <- factor(paste0(samples$plantline, "-", samples$time_point))
#' Set reference level
samples$condition <- relevel(samples$condition, "t89-0w")
#' Update sample info in the file
write.csv(samples, file= "~/Git/UPSCb/projects/aba/doc/sampleInfo_allSamples_BioRep.csv")

#' Create the dds object  
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = samples,
  design = ~ condition)

save(dds,file = "../DESeq2/DESeq2_object_t89aba.rda")

#' Check the size factors (i.e. the sequencing library size effect)
#' There is no big variation, a Variance Stabilizing Transformation can
#' be used (over a Relative Log Transformation)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count.table)
sizes
boxplot(sizes, main="Sequencing libraries size factor")

#' Perform the VST
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)
vst <- vst - min(vst)

save(vst,file = "library-size-normalized_variance-stabilized_data_t89aba.rda")
write.csv(vst,"library-size-normalized_variance-stabilized_data_t89aba.csv")

#' Validate the VST 
#' 
#' Visualize the corrected mean - sd
#' relationship. It is fairly linear,
#' meaning we can assume homoscedasticity.
#' The slight initial trend / bump is
#' due to genes having few counts in
#' a few subset of the samples and hence 
#' having a higher variability. This is
#' expected.
meanSdPlot(vst[rowSums(count.table)>0,])
any(rowSums(count.table) == 0)

#' # QC on the normalised data
#' 
#' ## PCA
#' 
#' First perform a Principal Component
#' Analysis (PCA) of the data
#'  to do a quick quality assessment; 
#' i.e. replicate should cluster
#' and the first 2-3 dimensions should 
#' be explainable by biological means.
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Plot the PCA 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(factor(samples$time_point))],
              pch=19)
legend("topright",pch=19,
       col=pal[1:3],
       cex = 0.8,
       legend=levels(factor(samples$time_point)))
par(mar=mar)

#' Then the first two dimensions
#' The first dimension separates the samples with the time points which is good
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(samples$time_point))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("topleft",
       pch=19,
       col=pal[1:3],
       cex = 0.9,
        legend=levels(factor(samples$time_point)))

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[c(5,6)][as.integer(factor(samples$plantline))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("topleft",
       pch=19,
       col=pal[c(5,6)],
       cex = 0.9,
       legend=levels(factor(samples$plantline)))

#' And the 2nd and 3rd dims
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(samples$time_point))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("top",
       pch=19,
       col=pal[1:3],
       cex = 0.9,
       legend=levels(factor(samples$time_point)))

plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[c(5,6)][as.integer(factor(samples$plantline))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("top",
       pch=19,
       col=pal[c(5,6)],
       cex = 0.9,
       legend=levels(factor(samples$plantline)))

#' There seems to be no outlier.  
#' 

#' ## Heatmap
#' The 1000 most variable genes are selected and plotted as a heatmap
dimnames(vst) <- dimnames(count.table)
sel <- order(apply(vst,1,sd),decreasing=TRUE)[1:1000]

heatmap.2(vst[sel,],labRow = NA,trace = "none", cexCol = 0.8 )
