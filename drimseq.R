library("tximport")
library("readr")
library("GenomicFeatures")
library("DRIMSeq")
library("txdbmaker")
library("BiocParallel")
library('ggplot2')

###trying BiocParallel() for time savings and adding threads_number

csv_dir <- "results"
samps <- read.csv(file.path(csv_dir, "script/samples.csv"))
head(samps)
samps$condition <- factor(samps$condition, levels=c("Control","Aza","Gua"))
quant_dir <- "results/salmon"
files <- file.path(quant_dir, samps$sample_id, "quant.sf")
names(files) <- samps$sample_id
head(files)

txi <- tximport(files, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")#tximport DTU checken nach countsFromAbundace = "dtuScaledTPM" statt scaled.TPM
cts <- txi$counts
cts <- cts[rowSums(cts)>0,]
rownames(cts) <- sub("\\|.*$", "", rownames(cts)) #fixing rownames to match with gtf file
head(cts)

gtf <- "data/ref/gencode.v50.chr_patch_hapl_scaff.annotation.gtf"
txdb.filename <- "gencode.50.annotation.sqlite"
if (file.exists(txdb.filename)){
    new_txdb <- loadDb(txdb.filename)
    txdf <- select(new_txdb, keys(new_txdb, "GENEID"),"TXNAME","GENEID")
    tab <- table(txdf$GENEID)
    txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
} else {
    txdb <- makeTxDbFromGFF(gtf)
    saveDb(txdb, txdb.filename)
}

sum(rownames(cts) %in% txdf$TXNAME)
sum(!rownames(cts) %in% txdf$TXNAME)
head(setdiff(rownames(cts), txdf$TXNAME))


#all(rownames(cts) %in% txdf$TXNAME)
stopifnot("Transkripte aus Salmon fehlen in der GTF" =
            all(rownames(cts) %in% txdf$TXNAME))
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
#all(rownames(cts) == txdf$TXNAME)
stopifnot(identical(rownames(cts), txdf$TXNAME))
counts <- data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts)

d <- dmDSdata(counts=counts, samples=samps)
methods(class=class(d))
counts(d[1,])[,1:4]
n_samps <- 9
n_samps.small <- 3
d <- dmFilter(d, min_samps_gene_expr=n_samps, min_gene_expr=10, min_samps_feature_expr=n_samps.small,min_feature_expr=10,min_samps_feature_prop=n_samps.small,min_feature_prop=0.1)#min_feature_expr pretty strict (trenscript does not make up more than 10% of gene's expression for at least n.small samlpes): else omit or lower min_feature_prop::swimming downstream
d
table(table(counts(d)$gene_id))
#DRIMSeq::samples(d)$condition <- factor(DRIMSeq::samples(d)$condition, levels = c("Control", "Aza", "Gua"))#should changes intercept/ref to Control
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)
design_full

n_workers <- as.integer(16)
bp <- MulticoreParam(workers = n_workers, RNGseed = 1)
set.seed(1)
d <- dmPrecision(d, design=design_full, BPPARAM=bp)
d <- dmFit(d, design=design_full, BPPARAM=bp)
d_aza <- dmTest(d, coef="conditionAza", BPPARAM=bp)#what to choose for 3 conditions? --> coef=c("ConditionAza","ConditionGua") als globaler Test, d_aza = dmTest(d,coef="ConditionAza") danach "ConditionGua" als einzelne Tests
#und ohne coef mit contrast=c(0,1,-1) zum Vergleich der Koeff von Aza vs Gua
d_gua <- dmTest(d, coef="conditionGua", BPPARAM=bp)
d_any <- dmTest(d, coef=c("conditionAza", "conditionGua"), BPPARAM=bp)
d_coef <- dmTest(d, contrast=c(0,1,-1), BPPARAM=bp)
res_aza <- DRIMSeq::results(d_aza)
head(res_aza)
res_aza.txp <- DRIMSeq::results(d_aza, level="feature")
head(res_aza.txp)

results <- list(res_gua = DRIMSeq::results(d_aza), res_gua_txp = DRIMSeq::results(d_aza, level="feature"), res_aza=DRIMSeq::results(d_gua), res_aza_txp=DRIMSeq::results(d_gua, level="feature"), res_coef=DRIMSeq::results(d_coef), res_coef_txp=DRIMSeq::results(d_coef, level="feature"), res_any=DRIMSeq::results(d_any))
no.na <- function(x) ifelse(is.na(x), 1, x)
out_dir <- "results/drimseq"
for (nam in names(results)) {
    res <- results[[nam]]
    write.csv(res, file.path(out_dir, paste0(nam, ".csv")), row.names = FALSE)
    res$pvalue <- no.na(res$pvalue)
    res <- res[order(res$adj_pvalue, res$pvalue), ]
    idx <- which(res$adj_pvalue < 0.05)[1]
    res[idx,]
    plot <- plotProportions(d, res$gene_id[idx], "condition")
    ggplot2::ggsave(file.path(out_dir, paste0(nam, "_proportions.png")),
                    plot = plot, width=8, height=5, dpi=300)
}
