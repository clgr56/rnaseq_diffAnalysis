library("DEXSeq")
library("DRIMSeq")
library("GenomicFeatures")
library("tximport")

gtf <- "data/ref/gencode.v50.chr_patch_hapl_scaff.annotation.gtf"

d <- readRDS('results/salmon/drimseqfilter_obj.rds', refhook=NULL)
sample.data <- DRIMSeq::samples(d)
####chat verbesserung!!!!!!!!!
# Sample-IDs als Row Names setzen
rownames(sample.data) <- sample.data$sample_id

# Sample-Reihenfolge an Count-Matrix anpassen
sample.data <- sample.data[
    colnames(count.data),
    ,
    drop = FALSE
]

# Condition als Faktor definieren
sample.data$condition <- factor(
    sample.data$condition,
    levels = c("Control", "Aza", "Gua")
)

# Kontrolle
stopifnot(
    identical(rownames(sample.data), colnames(count.data))
)
#!!!!!!!!!!

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
dxr.t = as.data.frame(dxr)

#columns <- c("featureID","groupID","pvalue") würde padj entfernen
dxr <- as.data.frame(dxr)
head(dxr)
dim(dxr.g[dxr.g$qval < 0.05,])
dim(dxr[!is.na(dxr$padj) & dxr$padj < 0.05, ])
#which(!mcols(dxd)$fullBetaConv) um problematische Zeile anzuzeigen
###chat vorschläge!!!!!!!!!!!!!!!!!!!!
nrow(dxr)
sum(!is.na(dxr$padj))
sum(dxr$padj < 0.05, na.rm = TRUE)

# Anteil signifikanter Transkripte
mean(dxr$padj < 0.05, na.rm = TRUE)

sig.transcripts <- subset(
    dxr,
    !is.na(padj) & padj < 0.05
)

# Anzahl signifikanter Transkripte pro Gen
sig.per.gene <- sort(
    table(sig.transcripts$groupID),
    decreasing = TRUE
)

head(sig.per.gene, 20)

sample.data
colnames(count.data)
rownames(sample.data)

# Prüfen, ob die Reihenfolge übereinstimmt
all(colnames(count.data) == rownames(sample.data))

# Anzahl eindeutiger Gene mit signifikanten Transkripten
length(unique(sig.transcripts$groupID))

# Anzahl signifikanter Transkripte pro Gen
sort(
    table(sig.transcripts$groupID),
    decreasing = TRUE
)[1:20]

# Anteil signifikanter Features
mean(dxr$padj < 0.05, na.rm = TRUE)
#!!!!!!!!!!!!!!!!!!!!!!!