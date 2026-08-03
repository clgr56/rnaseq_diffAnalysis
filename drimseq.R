library("tximport")
library("readr")
library("GenomicFeatures")
library("DRIMSeq")
library("txdbmaker")

csv_dir <- "results/"
samps <- read.csv(file.path(csv_dir, "samples.csv"))
head(samps)
quant_dir <- "results/salmon/"
files <- file.path(quant_dir, samps$sample_id, "quant.sf")
names(files) <- samps$sample_id
head(files)

txi <- tximport(files, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts)>0,]
head(cts)

gtf <- "ref/gencode.v50.chr_patch_hapl_scaff.annotation.gtf"
txdb.filename <- "gencode.50.annotation.sqlite"
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)
#load txdb
new_txdb <- loadDb(txdb.filename)
txdf <- select(txdb, keys(txdb, "GENEID"),"TXNAME","GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
