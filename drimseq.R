library("tximport")
library("readr")
library("GenomicFeatures")
library("DRIMSeq")
library("txdbmaker")

csv_dir <- "results"
samps <- read.csv(file.path(csv_dir, "script/samples.csv"))
head(samps)
quant_dir <- "results/salmon"
files <- file.path(quant_dir, samps$sample_id, "quant.sf")
names(files) <- samps$sample_id
head(files)

txi <- tximport(files, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts)>0,]
rownames(cts) <- sub("\\|.*$", "", rownames(cts)) #fixing rownames to match with gtf file
head(cts)
sum(rownames(cts) %in% txdf$TXNAME)
sum(!rownames(cts) %in% txdf$TXNAME)
head(setdiff(rownames(cts), txdf$TXNAME))


gtf <- "data/ref/gencode.v50.chr_patch_hapl_scaff.annotation.gtf"
txdb.filename <- "gencode.50.annotation.sqlite"
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)
#load txdb
new_txdb <- loadDb(txdb.filename)
txdf <- select(txdb, keys(txdb, "GENEID"),"TXNAME","GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]


all(rownames(cts) %in% txdf$TXNAME)
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)
counts <- data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts)

d <- dmDSdata(counts=counts, samples=samps)
methods(class=class(d))
counts(d[1,])[,1:4]
n_samps <- 9
n_samps.small <- 3
d <- dmFilter(d, min_samps_gene_expr=n, min_gene_expr=10, min_samps_feature_expr=n.small,min_feature_expr=10,min_samps_feature_prop=n.small,min_feature_prop=0.01)#min_feature_expr pretty strict (trenscript does not make up more than 10% of gene's expression for at least n.small samlpes): else omit or lower min_feature_prop::swimming downstream
d
table(table(counts(d)$gene_id))
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)
design_full
#set.seed(1)
#d <- dmPrecision(d, design=design_full)
#d <- dmFit(d, design=design_full)
#d <- dmTest(d, coef="")#what to choose for 3 conditions
