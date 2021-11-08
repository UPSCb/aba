#!/mnt/aspnas/sw/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
    stop('usage: runDESeq2.R <sample file> <htseq dir> <output dir>')
}

library(DESeq2)
source('~niklasm/git/UPSCb/src/R/plotDispLSD.R')

sample_file = args[1]
htseq_dir = args[2]
output_dir = args[3]

sample = strsplit(tail(strsplit(sample_file, '/')[[1]], n=1), '\\.')[[1]][1]

dir.create(file.path(output_dir, sample))
output_prefix = paste(output_dir, sample, sample, sep='/')

sampleTable = read.table(sample_file, header=TRUE)
condition_names = colnames(sampleTable)[3:4]
#conditions = levels(sampleTable[,3])
conditions = c(as.character(sampleTable[1,3]), as.character(sampleTable[nrow(sampleTable),3]))

if (length(levels(sampleTable[,3])) > 2) {
    stop('you can compare maximum 2 conditions')
} else {
    cat(paste('Comparing', conditions[1], 'against', conditions[2], "\n", sep=' '))
}

ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable, htseq_dir, ~ condition)
# Make sure that the log ratios are calculated as intended
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, levels=c(conditions[2], conditions[1]))

dds = DESeq(ddsHTSeq)

res = results(dds)
res = res[order(res$padj),]

png(paste(output_prefix, '_dispersions.png', sep=''), width=800, height=600)
plotDispEsts(dds)
dev.off()

png(paste(output_prefix, '_MA.png', sep=''), width=800, height=600)
plotMA(dds)
dev.off()

# Write the DEGs to file
write.table(as.data.frame(res), file=paste(output_prefix, '_de.txt', sep=''),
    sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)

library(RColorBrewer)
library(gplots)
library(genefilter)

dds.norm = counts(dds, normalized=TRUE)
dists = dist(t(dds.norm))
mat = as.matrix(dists)

hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

png(paste(output_prefix, "_sampledist.png", sep=""), width=800, height=600)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))
dev.off()

select = order(rowVars(dds.norm), decreasing=TRUE)[seq_len(500)]
fac = factor(apply(as.data.frame(colData(dds))[, 'condition', drop=FALSE], 1, paste))
colours = c('olivedrab', 'cornflowerblue')

pca = prcomp(t(dds.norm[select,]))

pca_perc = round(summary(pca)$importance[2,]*100)

png(paste(output_prefix, "_samplepca.png", sep=""), width=800, height=600)
xyplot(PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=2, col=colours, aspect="iso",
    main=draw.key(key=list(rect=list(col=colours), text=list(levels(fac)), rep=FALSE)),
    xlab=paste('PC1 (', pca_perc[1], '%)', sep=''), ylab=paste('PC2 (', pca_perc[2], '%)', sep=''))
dev.off()

# Write a description to file
write.table(as.data.frame(mcols(res, use.names=TRUE)), file=paste(output_prefix, '_description.txt', sep=''),
    sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
