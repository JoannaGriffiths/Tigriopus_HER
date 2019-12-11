source('http://bioconductor.org/biocLite.R')
biocLite('DESeq2')
biocLite("edgeR")
biocLite("limma")
biocLite("statmod")
biocLite("affycoretools")
biocLite("ReportingTools")
library("edgeR")
library("statmod")
library('DESeq2')

setwd("~/Documents/Tigriopus_whole_genome_seq/gene_exp/DeSeq2")
countdata_1 <- read.delim("Tig_RSEM_merged_matrix", header=TRUE, row.names="Library") #read in the file of the count data and call it countdata, row.names tells the name in the top left cell, the gene names

countdata_1 = round(countdata_1)

keep <- rowSums(cpm(countdata_1)>1) >=4
# keeps rows (genes) where at least 4 columns (libraries) have at least 1 count per million. This means that if a gene is only expressed in say one treatment (which has three replicates), this gene will not be thrown out of the analysis

countdata<- countdata_1[keep,] #formatting for organizing the kept rows that summed to at least 1 cpm in the step above

countdata=round(countdata)

coldata <- read.delim("column_Tig_all.txt", header=TRUE, row.names=1)

#select one from below:
countdata = countdata[, c(1,2,3,6,7,10)] #select all heat shocked
countdata = countdata[, c(4,5,8,9,11,12)] #select all controls
countdata = countdata[, c(1,2,3,4,5,11)] #select all unselected lines
countdata = countdata[, c(6,7,8,9,10,12)] #select all selected lines


ddsFullCountTable <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design = ~ Selected)

ddsFull <- DESeq(ddsFullCountTable) # this is the analysis!
head(ddsFull)
res=results(ddsFull)
res
resOrdered=res[order(res$padj),]
head(resOrdered)
sum(res$padj<0.05, na.rm=TRUE)

write.table(resOrdered, "Tig_output_DEG_all_sel_vs_unsel.txt")

rld <-rlogTransformation(ddsFull)

write.table(assay(rld), "DEG_logtransformed_Tig_all_sel_vs_unsel", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

head(assay(rld))

hist(assay(rld))

#install.packages("gplots")
library(gplots) 
#install.packages("RColorBrewer")
library("RColorBrewer") 
library( "genefilter" )


topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 40)


quartz()
heatmap.2(assay(rld)[topVarGenes, ], scale="row",
          trace="none", dendrogram="both", key=TRUE, keysize = 1.5, margins =c(3,11), density.info = "density",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = heatmap.2( assay(rld)[ topVarGenes, ], # category labels
                           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),  # color key
                           lty= 1,             # line style
                           lwd = 10))          # line width




