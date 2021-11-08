#' ---
#' title: "ABA T89 - Biological QA (Kallisto)"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---

#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/aspseq/rbhalerao/ABA/Potra01/")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/aspseq/rbhalerao/ABA/Potra01/")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(stringr))


#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Create palettes
pal <- brewer.pal(8,"Dark2")
pal12 <- brewer.pal(12,"Paired")
tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' Register the default plot margin
mar <- par("mar")

#' # Raw data
#' ## Loading  
#' ### Original data
orig <- list.files("kallisto", 
                    recursive = TRUE, 
                    pattern = "abundance.tsv",
                    full.names = TRUE)

#' name them
names(orig) <- sub("_sortmerna_trimmomatic*","",
                    sapply(strsplit(orig, "/"), .subset, 2))

#' Create sample information
samples <- data.frame(file_name = str_extract(orig, "[0-9]+w.+001"),
                      type = rep("T89", 11),
                      time = str_extract(orig, "[0-9]+w"),
                      replicate = gsub("_", "", str_extract(orig, "_[0-9]_")))

#' Reorder the sample data.frame according to the way
#' the results were read in and accounting for tech reps
samples <- samples[match(names(orig),samples$file_name),]

#' Read the tx2gene translation
tx2gene <- read.table("/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/kallisto/tx2gene.tsv",
                      header=FALSE)

#' Read the expression at the transcript level
tx <- suppressMessages(tximport(files = orig, 
                                type = "kallisto", 
                                txOut = TRUE))
kt <- round(tx$counts)

kg <- round(summarizeToGene(txi=tx,tx2gene=tx2gene)$counts)

#' Some genes from Kallisto index are not present in the tx2gene. Update tx2gene and repeat the analysis!  
#' ## Raw Data QC analysis  
#' ### Original data  

#' Check how many genes are never expressed
sel <- rowSums(kg) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(kg),digits=1),
        sum(sel),
        nrow(kg))

#' The cumulative gene coverage is as expected
plot(density(log10(rowMeans(kg))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by time. 
plot.multidensity(lapply(1:ncol(kg),function(k){log10(kg)[,k]}),
                  col=c(1,pal)[as.integer(samples$time)],
                  legend.x="topright",
                  legend=levels(samples$time),
                  legend.col=c(1,pal)[1:nlevels(samples$time)],
                  #legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' Most samples show similar trend. Some are shifted to the left, probably they were not sequenced so deeply.  

#' ## Raw data export
dir.create(file.path("analysis","kallisto"),showWarnings=FALSE, recursive = TRUE)
write.csv(kg,file="analysis/kallisto/T89_raw-unnormalised-gene-expression_data.csv")

#' # Data normalisation 
#' ## Preparation  
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate  

#' ### Original data  
#' Before creating dds object, make sure the order of the factors is correct (0W set as the first one) for later
#' DE analysis.  
samples$time <- relevel(samples$time,"0w")

#' Create the dds object, without giving any prior on the design
dds.kg <- DESeqDataSetFromMatrix(
  countData = kg,
  colData = samples,
  design = ~time)

#' Check the size factors (i.e. the sequencing library size effect)
dds.kg <- estimateSizeFactors(dds.kg)
sizes.kg <- sizeFactors(dds.kg)
names(sizes.kg) <- colnames(kg)
pander(sizes.kg)
boxplot(sizes.kg, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
#' At the gene level
vsd.kg <- varianceStabilizingTransformation(dds.kg, blind=TRUE)
vst.kg <- assay(vsd.kg)
vst.kg <- vst.kg - min(vst.kg)

#' Validate the VST 
meanSdPlot(vst.kg[rowSums(kg)>0,])

#' Export the vst
write.csv(vst.kg,"analysis/kallisto/T89_library-size-normalized_variance-stabilized_gene-expression_data.csv")

#' # QC on the normalised data
#' 
#' ## PCA
".pca" <- function(vst,fact,lgd="bottomright",pal=brewer.pal(8,"Dark2")){
  pc <- prcomp(t(vst))
  
  percent <- round(summary(pc)$importance[2,]*100)
  
  #' ### 3 first dimensions
  mar=c(5.1,4.1,4.1,2.1)
  scatterplot3d(pc$x[,1],
                pc$x[,2],
                pc$x[,3],
                xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
                ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
                zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
                color=pal[as.integer(fact)],
                pch=19)
  legend(lgd,pch=19,
         col=pal[1:nlevels(fact)],
         legend=levels(fact))
  par(mar=mar)
}

#' ### Time
.pca(vst.kg,factor(samples$time),lgd="topright",pal=c(1,pal))

#' ### 1st and 2nd dims
pc <- prcomp(t(vst.kg))
percent <- round(summary(pc)$importance[2,]*100)

#' Time

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal12[as.integer(factor(samples$time))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("bottomleft",bty="n",col=pal12,levels(samples$time),pch=19)

#' ### 2nd and 3rd dims

#' Time
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal12[as.integer(factor(samples$time))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topleft",bty="n",col=pal12,levels(samples$time),pch=19)

#' Replicate
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=c(1,pal)[as.integer(samples$time)],
     pch=c(1:nlevels(samples$replicate)),
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("bottomleft",bty="n",col=c(1,pal), levels(samples$time), pch = 19)
legend("topleft", bty = "n", pch=c(1:nlevels(samples$replicate)), levels(samples$replicate))

#' ## Heatmap
sel <- featureSelect(vst.kg, samples$time,exp = 3,nrep = 2)
message(sprintf("We select %s genes for the heatmap",sum(sel)))
heatmap.2(vst.kg[sel,],trace = "none",labRow = FALSE, 
          col=hpal, ColSideColors = c(1,pal)[as.integer(samples$time)],
          scale="row")

legend("bottomright", bty="n", col=c(1,pal), levels(samples$time),pch=19)

#' ### Sample Hierarchical clustering
hc <- hclust(dist(t(vst.kg[sel,])))
plot(hc, labels=samples$time,
     main = "Hierarchical clustering",cex=0.7)

#' # Merge the technical replicates
counts <- do.call(
  cbind,
  lapply(split.data.frame(t(kg),
                          paste0(samples$time, "-", samples$replicate)),
         colSums))

csamples <- samples[,-1]
csamples <- csamples[match(colnames(counts),paste0(samples$time, "-", samples$replicate)),]
write.csv(csamples,file="~/Git/UPSCb/projects/aba/doc/biological-samples_T89.csv")

#' ## QC 
#' The cumulative gene coverage is as expected
plot(density(log10(rowMeans(counts))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by time. 
plot.multidensity(lapply(1:ncol(counts),function(k){log10(counts)[,k]}),
                  col=c(1,pal)[as.integer(csamples$time)],
                  legend.x="topright",
                  legend=levels(csamples$time),
                  legend.col=c(1,pal)[1:nlevels(csamples$time)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' Write the count table (tech. rep. combined)
write.csv(counts,file="analysis/kallisto/T89_raw-unnormalised-gene-expression-tech-rep-combined_data.csv")

#' # Data normalisation 
#' ## Preparation
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate  

#' ### Original data  
#' Create the dds object, without giving any prior on the design
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = csamples,
  design = ~time)

#' Check the size factors (i.e. the sequencing library size effect)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(counts)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")
pander(sort(sizes))

#' ## Variance Stabilising Transformation
#' At the gene level
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' Validate the VST 
meanSdPlot(vst[rowSums(vst)>0,])

#' Export the vst
write.csv(vst,"analysis/kallisto/T89_library-size-normalized_variance-stabilized_gene-expression-tech-rep-combined_data.csv")

#' # QC on the normalised data
#' 
#' ## PCA
".pca" <- function(vst,fact,lgd="topleft",pal=brewer.pal(8,"Dark2")){
  pc <- prcomp(t(vst))
  
  percent <- round(summary(pc)$importance[2,]*100)
  
  #' ### 3 first dimensions
  mar=c(5.1,4.1,4.1,2.1)
  scatterplot3d(pc$x[,1],
                pc$x[,2],
                pc$x[,3],
                xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
                ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
                zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
                color=pal[as.integer(fact)],
                pch=19)
  legend(lgd,pch=19,
         col=pal[1:nlevels(fact)],
         legend=levels(fact))
  par(mar=mar)
}

#' ### Time
.pca(vst,factor(csamples$time),lgd="topright",pal=c(1,pal))

#' ### 1st and 2nd dims
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Time
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal12[as.integer(factor(csamples$time))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("bottomleft",bty="n",col=pal12,levels(csamples$time),pch=19)

#' ### 2nd and 3rd dims  
#' Time
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal12[as.integer(factor(csamples$time))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("bottomleft",bty="n",col=pal12,levels(csamples$time),pch=19)

#' ## Heatmap
sel <- featureSelect(vst,csamples$time,exp = 3,nrep = 2)
heatmap.2(vst[sel,],trace = "none",labRow = FALSE, 
          col=hpal, ColSideColors = c(1,pal)[as.integer(csamples$time)],
          scale="row")

#' ### Sample Hierarchical clustering
hc <- hclust(dist(t(vst[sel,])))
plot(hc, labels=csamples$triplicate,
     main = "Hierarchical clustering",cex=0.7)

#'  
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
