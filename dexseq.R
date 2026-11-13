library("DEXSeq")
library("DRIMSeq")
library(GenomicFeatures)
library(tximport)

txi <- tximport(files, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")#tximport DTU checken nach countsFromAbundace = "dtuScaledTPM" statt scaled.TPM
cts <- txi$counts
cts <- cts[rowSums(cts)>0,]
rownames(cts) <- sub("\\|.*$", "", rownames(cts)) #fixing rownames to match with gtf file
head(cts)

gtf <- "data/ref/gencode.v50.chr_patch_hapl_scaff.annotation.gtf"
txdb.filename <- "gencode.50.annotation.sqlite"
new_txdb <- loadDb(txdb.filename)
txdf <- select(new_txdb, keys(new_txdb, "GENEID"),"TXNAME","GENEID")
head(txdf)
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

counts <- data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts)

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
