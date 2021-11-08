#' ---
#' title: "Differential Expression analysis of samples from T89 and abi1 plants in dormancy cycle"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---     
          
#' # Setup
#' # Environment
#' Set the working dir
setwd("/mnt/picea/projects/aspseq/rbhalerao/ABA/Potrx01/analysis/all_samples/DESeq2/")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/aspseq/rbhalerao/ABA/Potrx01/analysis/all_samples/DESeq2/")
#' ```

#' Libs
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(VennDiagram))

#' Helper files
source("~/Git/UPSCb/src/R/plotMA.R")
source("~/Git/UPSCb/src/R/volcanoPlot.R")

#' Load saved data
load("DESeq2_object_t89aba.rda")

#' ########### should I use the same value or the usual settings of lcf 0.5 and padj 0.01? How were these values selected?
#' set the Fold change cut off value (2017-03-16 originally it was 1 but change it now ) # 0.4 means about 32% change in expression 
RNAexpr_log2Fold_cutoff = 0.4

#' set the cutoff value
set.alpha = 1e-6

# set alpha value for MA plot
s.alpha = 1e-6

# create an object with the actual cut off values for saving that also into the file names for the selected genes in csv files
lfch.cutoff_padj <-  paste0("lfcoff",RNAexpr_log2Fold_cutoff, "_padjcoff_e6")
# Give names for samples, time points
names.comparison <- c("0W_abi1_wt","6WSD_abi1_wt", "10WSD_abi1_t89")

###########################################



#' Setup graphics
pal=brewer.pal(8,"Dark2")
mar <- par("mar")

#' # Process
#' run DESeq2 (on the condition - i.e. combining replicates)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

#' plot the disp. est.
plotDispEsts(dds)

################## why running a specific test?
#' #' run the test
#' dds <- nbinomWaldTest(dds)
#################  

#' Run DE
dds <- DESeq(dds)

#' Extract DEGs between WT and abi1 mutant in each of the time points
res_abi1WT_0 <- results(dds, contrast = c("condition", "abi1-0w", "t89-0w"), filter = rowMedians(counts(dds)))
res_abi1WT_6 <- results(dds, contrast = c("condition", "abi1-6w", "t89-6w"), filter = rowMedians(counts(dds)))
res_abi1WT_10 <- results(dds, contrast = c("condition", "abi1-10w", "t89-10w"), filter = rowMedians(counts(dds)))


#' plot the MA and volcano plot
MA.plot <- function(res, salpha){
  #par(mar=c(5.1,5.1,5.1,0.5))
  plotMA(res,alpha= salpha)
  #par(mar=mar)
}


# Create a list from the results  

res.list<-list(res_abi1WT_0, res_abi1WT_6, res_abi1WT_10)
names(res.list) <- c("0W_abi1_wt", "6W_abi1_wt", "10W_abi1_wt")

# generate MA plots
par(mar=c(5.1,5.1,5.1,0.5))
lapply(res.list, MA.plot, s.alpha)

#'create volcano plots with fixed x and y values and set points color
#'
volcano.XYc <- function (object, alpha, point.col){
  sel <- ! is.na(object$padj)
  sel2 <- object$padj[sel]<=alpha 
  sel3up <- object$log2FoldChange[sel2]>= RNAexpr_log2Fold_cutoff 
  sel3dw <- object$log2FoldChange[sel2]<= -RNAexpr_log2Fold_cutoff
  
  ## plot
  heatscatter(object$log2FoldChange[sel],
              -log10(object$padj[sel]),
              main="Volcano",
              xlab=elementMetadata(object)$description[2], 
              ylab="- log(10) adj. p-value",
              xlim = c(-12,12),
              ylim = c(0,300)
              #,add.contour=TRUE,colpal="spectral"
  )
  
  ## legend
  legend("topleft",bty="n",paste("cutoff @",alpha),lty=2,col="gray")
  
  # choose the color of the points
  if (point.col == "blue")
  { col.1 <- "lightblue"
    col.2  <- "dodgerblue3"
  } else if (point.col == "green"){
    col.1 <- "darkseagreen3"
    col.2  <- "darkseagreen4"
  }
  points(object$log2FoldChange[sel][sel2],-log10(object$padj[sel][sel2]),col=col.1 ,pch=19)
  points(object$log2FoldChange[sel][sel2],-log10(object$padj[sel][sel2]),col=col.2 ,pch=19,cex=0.5)
  abline(h=-log10(alpha),lty=2,col="gray")
}

#' Create volcano plot from all the results
par(mfrow=c(2,2))

lapply(res.list, volcano.XYc, s.alpha, "blue")

#' It's actually interesting that the adjusted p-value are so
#' extreme (log odds of 200).
#' If we think about the "q" value (i.e. the odds that we make an
#' observation) we might expect maybe 1% of the genes to be affected
#' and that means that an adj. p-value of 0.001 (50 times less than a 0.05 
#' value for a 50/50 % odd) might not still be sufficient to warrant a high
#' chance of a real effect. Even a 10-6 cutoff shows a lot of "significant"
#' genes on the volcanoplot.
#' 

#' Count the number of genes between the stages   
#' for that the function 
diff.gene.number<- function (re, valpha){
  alpha=valpha
  sel <- re$padj<=alpha
  sprintf("There are %s genes differentially expressed at a %s cutoff and %s 
genes are positively differentially expressed.  That is %.2f percentage of the changing ones in that comparision",
          sum(sel,na.rm=TRUE),alpha ,sum(sign(re$log2FoldChange[sel])==1,na.rm=TRUE)
          ,((sum(sign(re$log2FoldChange[sel])==1,na.rm=TRUE))/(sum(sel,na.rm=TRUE)))*100)
  }

#' set the cutoff value
set.alpha = 1e-2
#'and the number of differentially expressed genes :
lapply(res.list,diff.gene.number , set.alpha)

set.alpha = 1e-6
#'and the number of differentially expressed genes :
lapply(res.list,diff.gene.number , set.alpha)


# go with the gene list selections  

#function for gene selection with increased or decreased expression values 
RNA.expr.g.sel <- function(res, incr_or_decr) {

  if (incr_or_decr == "incr")
  {g.res <- rownames(res[res$padj<set.alpha & ! is.na(res$padj) & res$log2FoldChange > {RNAexpr_log2Fold_cutoff},])
  }
  else if (incr_or_decr == "decr")
  {g.res <- rownames(res[res$padj<set.alpha & ! is.na(res$padj) & res$log2FoldChange < -{RNAexpr_log2Fold_cutoff},])
  }
# transform the names to matrix and / take away the ".0 "
  g.res  <- sub("\\.0","",as.matrix(g.res))
  return (g.res)
}

# select genes with increased or dicreased expressions between plant types
res_decr <- lapply(res.list, RNA.expr.g.sel, "decr")
res_incr <- lapply(res.list, RNA.expr.g.sel, "incr")

# number of selected DEGs
elementNROWS(res_decr)
elementNROWS(res_incr)

#' save the gene lists into csv file ## hoe to get rid of "V1" (col.names = FALSE does not work, although it is supposed to)?
lapply(seq_along(names(res_decr)), function(x){
  write.csv(res_decr[[x]], file = paste0(names(res_decr)[[x]], "_", lfch.cutoff_padj, "_decr.csv"), row.names = FALSE, quote = FALSE)
})

lapply(seq_along(names(res_incr)), function(x){
  write.csv(res_incr[[x]], file = paste0(names(res_incr)[[x]], "_", lfch.cutoff_padj, "_incr.csv"), row.names = FALSE, quote = FALSE)
})

#' # Session Info
sessionInfo()