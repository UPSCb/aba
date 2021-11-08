#' ---
#' title: "ABA experiment sample comparison"
#' author: "Nicolas Delhomme & Pal Miskolczi"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/aspseq/ABA/Potrx01")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/aspseq/ABA/Potrx01")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(vsn))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Create a palette
pal <- brewer.pal(8,"Dark2")

#' Register the default plot margin
mar <- par("mar")

#' Read the sample information
  samples <- read.csv("~/Git/UPSCb/projects/aba/doc/samples.csv")

#' Read the HTSeq files in a matrix
res <- mclapply(dir("htseq",pattern=".*t89.*\\.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=max(mcaffinity()))

names(res) <- sub("_sortmerna.*\\.txt","",dir("htseq",pattern=".*t89.*\\.txt"))

#' Reorder the sample data.frame according to the way
#' the results were read in
samples <- samples[match(names(res),samples$SampleName),]

#' Raw Data QC analysis
addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]
dir.create(file.path("analysis","HTSeq"),recursive = TRUE, showWarnings = FALSE)
write.csv(count.table,"analysis/HTSeq/raw-unormalised-data.csv")

#' Extract the HTSeq stat lines
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

#' Convert them into percentages
apply(count.stats,2,function(co){round(co*100/sum(co))})

#' Plot the stats
#' 
#' There are no outliers
col <- pal[1:nrow(count.stats)]
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=col,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,4e+6),cex.names=.6)
legend("top",fill=col,legend=gsub("_"," ",rownames(count.stats)),bty="n",cex=0.8)
par(mar=mar)

#' The average percentage of aligned reads is 75%
boxplot(unlist(count.stats["aligned",]/colSums(count.stats)),
        main="aligned reads",ylab="percent aligned",ylim=c(0,1))

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
#' The cumulative coverage is as expected, around ???X
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by time points.
#' 
#' The observed distribution is relatively similar, a few samples show slight
#' diffences and there seem to be a trend of more average expression for
#' the  compared to the.
plot.multidensity(log10(count.table),
                  col=pal[as.integer(factor(samples$time_point))],
                  legend.x="topright",
                  legend=c("0W","6W","10W"),
                  legend.col=pal[1:3],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' # Data normalisation 
#'  For visualization, the data is
#' submitted to a variance stabilization
#' transformation using DESeq2. The 
#' dispersion is estimated independently
#' of the sample tissue and replicate

#' Create the dds object
conditions <- colnames(count.table)
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions),
  design = ~ condition)

#' Check the size factors (i.e. the sequencing library size effect)
#' There is no big variation, a Variance Stabilizing Transformation can
#' be used (over a Relative Log Transformation)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count.table)
sizes
boxplot(sizes, main="Sequencing libraries size factor")

#' Perform the VST
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)
vst <- vst - min(vst)
write.csv(vst,"analysis/HTSeq/library-size-normalized_variance-stabilized_data.csv")

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
meanSdPlot(vst[rowSums(count.table)>0,], ylim = c(0,2.5))

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
#' There seem to be 2 outliers, one ZE and one FMG
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(factor(samples$time_point))],
              pch=19)
legend("topleft",pch=19,
       col=pal[1:3],
       legend=levels(factor(samples$time_point)))
par(mar=mar)

#' Then the first two dimensions
#' The first dimension separates the 2 tissue which is good,
#' whereas the 2nd dimension seem to separate the outlier from
#' the rest
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(samples$time_point))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

#' And the 2nd and 3rd dims
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(samples$time_point))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

sprintf("The sample %s seem to be an outlier", samples$Seq_laneID[pc$x[,2] > 50])

#' ## Heatmap
#' The 1000 most variable genes are selected and plotted as a heatmap
dimnames(vst) <- dimnames(count.table)
sel <- order(apply(vst,1,sd),decreasing=TRUE)[1:1000]

#' First with their SciLife ID

heatmap.2(vst[sel,],labRow = NA,trace = "none",cexCol = 0.6 )   

#' The distinction between time points  are very clear in these 1000 first genes
#' 
#' 
#' # Conclusion
#' The raw quality of the data appears very good. The data normalisation
#' gives satisfying results (as far as the VST fit is concerned). The PCA and the
#' heatmap identifies no real outliers. 
