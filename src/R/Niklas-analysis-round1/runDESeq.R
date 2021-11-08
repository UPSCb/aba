#!/mnt/aspnas/sw/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
    stop('usage: runDESeq.R <sample file> <htseq dir> <output dir>')
}

library(DESeq)
source('~niklasm/git/UPSCb/src/R/plotDispLSD.R')

sample_file = args[1]
htseq_dir = args[2]
output_dir = args[3]

sample = strsplit(tail(strsplit(sample_file, '/')[[1]], n=1), '\\.')[[1]][1]

dir.create(file.path(output_dir, sample))
output_prefix = paste(output_dir, sample, sample, sep='/')

sampleTable = read.table(sample_file, header=TRUE)
condition_names = colnames(sampleTable)[3:4]
conditions = levels(sampleTable[,3])

if (length(conditions) > 2) {
    stop('you can compare maximum 2 conditions')
} else {
    cat(paste('Comparing', conditions[1], 'against', conditions[2], "\n", sep=' '))
}

cds = newCountDataSetFromHTSeqCount(sampleTable[, 1:3], htseq_dir)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

png(paste(output_prefix, '_dispersions.png', sep=''), width=800, height=600)
plotDispLSD(cds)
dev.off()

res = nbinomTest(cds, conditions[2], conditions[1])

png(paste(output_prefix, '_MA.png', sep=''), width=800, height=600)
plotMALSD(res)
dev.off()

# Make the column names more clear, what is A/B?
decols = colnames(res)
decols[3] = paste("baseMean", conditions[2], sep='.')
decols[4] = paste("baseMean", conditions[1], sep='.')
decols[5] = paste("foldChange(", conditions[1], "/", conditions[2], ")", sep="")

# Write the DEGs to file
write.table(res[order(res$padj), ], file=paste(output_prefix, '_de.txt', sep=''),
    sep='\t', quote=FALSE, row.names=FALSE, col.names=decols)

# cdsFull = newCountDataSetFromHTSeqCount(sampleTable, htseq_dir)
# cdsFull = estimateSizeFactors(cdsFull)
# cdsFull = estimateDispersions(cdsFull)

# cdsFullBlind = estimateDispersions(cdsFull, method="blind")
# vsdFull = varianceStabilizingTransformation(cdsFullBlind)

# library(RColorBrewer)
# library(gplots)

# hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

# dists = dist(t(exprs(vsdFull)))
# mat = as.matrix(dists)
# rownames(mat) = colnames(mat) = rownames(pData(cdsFullBlind))

# # Sample heatmap
# png(paste(output_prefix, "_sampledist.png", sep=''), width=800, height=600)
# heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))
# dev.off()

# # Sample PCA
# png(paste(output_prefix, "_samplepca.png", sep=''), width=800, height=600)
# plotPCA(vsdFull, intgroup=condition_names)
# dev.off()

library(RColorBrewer)
library(gplots)
library(genefilter)

cds.norm = counts(cds, normalized=TRUE)
dists = dist(t(cds.norm))
mat = as.matrix(dists)

hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

png(paste(output_prefix, "_sampledist.png", sep=""), width=800, height=600)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))
dev.off()

select = order(rowVars(cds.norm), decreasing=TRUE)[seq_len(500)]
fac = factor(apply(pData(cds)[, 'condition', drop=FALSE], 1, paste))
colours = c('olivedrab', 'cornflowerblue')

pca = prcomp(t(cds.norm[select,]))

pca_perc = round(summary(pca)$importance[2,]*100)

png(paste(output_prefix, "_samplepca.png", sep=""), width=800, height=600)
xyplot(PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=2, col=colours, aspect="iso",
    main=draw.key(key=list(rect=list(col=colours), text=list(levels(fac)), rep=FALSE)),
    xlab=paste('PC1 (', pca_perc[1], '%)', sep=''), ylab=paste('PC2 (', pca_perc[2], '%)', sep=''))
dev.off()
