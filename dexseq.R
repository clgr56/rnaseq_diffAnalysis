library("DEXSeq")
library("DRIMSeq")
library(GenomicFeatures)
library(tximport)

gtf <- "data/ref/gencode.v50.chr_patch_hapl_scaff.annotation.gtf"

d <- readRDS('results/salmon/drimseqfilter_obj.rds', refhook=NULL)
sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))

dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + condition:exon,
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet=TRUE)
dxd <- testForDEU(dxd, reducedModel=~sample + exon)

dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)
dxr.t = as.data.frame(dxr[, c("featureID","groupID","pvalue")])

columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])
head(dxr)
dim(dxr.g[dxr.g$qval < 0.05,])
dim(dxr[dxr$padj < 0.05,])
