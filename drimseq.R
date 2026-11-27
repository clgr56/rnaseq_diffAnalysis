library("tximport")
library("readr")
library("GenomicFeatures")
library("DRIMSeq")
library("txdbmaker")
library("BiocParallel")
library("ggplot2")

###trying BiocParallel() for time savings and adding threads_number

csv_dir <- "script"
samps <- read.csv(file.path(csv_dir, "samples.csv"))
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
    txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
    tab <- table(txdf$GENEID)
    txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
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
saveRDS(d, 'results/salmon/drimseqfilter_obj.rds', compress=TRUE)

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

# results <- list(res_aza = DRIMSeq::results(d_aza),
#                 res_aza_txp = DRIMSeq::results(d_aza, level="feature"),
#                 res_gua=DRIMSeq::results(d_gua), 
#                 res_gua_txp=DRIMSeq::results(d_gua, level="feature"), 
#                 res_coef=DRIMSeq::results(d_coef), 
#                 res_coef_txp=DRIMSeq::results(d_coef, level="feature"), 
#                 res_any=DRIMSeq::results(d_any))
# no.na <- function(x) ifelse(is.na(x), 1, x)
# out_dir <- "results/drimseq"
# for (nam in names(results)) {
#     res <- results[[nam]]
#     write.csv(res, file.path(out_dir, paste0(nam, ".csv")), row.names = FALSE)
#     res$pvalue <- no.na(res$pvalue)
#     res <- res[order(res$adj_pvalue, res$pvalue), ]
#     idx <- which(res$adj_pvalue < 0.05)[1]
#     res[idx,]
#     plot <- plotProportions(d, res$gene_id[idx], "condition")
#     ggplot2::ggsave(file.path(out_dir, paste0(nam, "_proportions.png")),
#                     plot = plot, width=8, height=5, dpi=300)
# }


# ==================================================
# 1. Collect DRIMSeq results
# ==================================================

results <- list(
    Aza_vs_Control = list(
        gene = DRIMSeq::results(d_aza),
        feature = DRIMSeq::results(d_aza, level = "feature")
    ),
    Gua_vs_Control = list(
        gene = DRIMSeq::results(d_gua),
        feature = DRIMSeq::results(d_gua, level = "feature")
    ),
    Aza_vs_Gua = list(
        gene = DRIMSeq::results(d_coef),
        feature = DRIMSeq::results(d_coef, level = "feature")
    )
)

# ==================================================
# 2. Output directory
# ==================================================

out_dir <- "results/drimseq"
dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ==================================================
# 3. Function: save results and plot top 10 DTU
# ==================================================

plot_top10_dtu <- function(
    result_list,
    comparison,
    d,
    top_n = 10,
    padj_cutoff = 0.05
) {

    cat("\n====================================\n")
    cat("Comparison:", comparison, "\n")
    cat("====================================\n")

    gene_res <- result_list$gene
    feature_res <- result_list$feature

    # --------------------------------------------------
    # Save complete results
    # --------------------------------------------------

    write.csv(
        gene_res,
        file.path(
            out_dir,
            paste0(comparison, "_gene_results.csv")
        ),
        row.names = FALSE
    )

    write.csv(
        feature_res,
        file.path(
            out_dir,
            paste0(comparison, "_transcript_results.csv")
        ),
        row.names = FALSE
    )

    # --------------------------------------------------
    # Select significant transcript features
    # --------------------------------------------------

    feature_res <- feature_res[
        !is.na(feature_res$adj_pvalue) &
        feature_res$adj_pvalue < padj_cutoff,
        ,
        drop = FALSE
    ]

    feature_res <- feature_res[
        order(feature_res$adj_pvalue, feature_res$pvalue),
        ,
        drop = FALSE
    ]

    # Select top N
    top_features <- head(feature_res, top_n)

    if (nrow(top_features) == 0) {
        cat("No significant DTU features found.\n")
        return(NULL)
    }

    cat(
        "Significant transcript features:",
        nrow(feature_res),
        "\n"
    )

    cat(
        "Plotting top",
        nrow(top_features),
        "transcript features.\n"
    )

    # --------------------------------------------------
    # Save top 10 transcript results
    # --------------------------------------------------

    write.csv(
        top_features,
        file.path(
            out_dir,
            paste0(comparison, "_top", top_n, "_DTU.csv")
        ),
        row.names = FALSE
    )

    # --------------------------------------------------
    # Identify unique genes containing top transcripts
    # --------------------------------------------------

    top_genes <- unique(top_features$gene_id)

    # --------------------------------------------------
    # Plot proportions for each gene
    # --------------------------------------------------

    for (gene in top_genes) {

        # Get transcript counts for this gene
        gene_counts <- DRIMSeq::counts(d)

        gene_counts <- gene_counts[
            gene_counts$gene_id == gene,
            ,
            drop = FALSE
        ]

        # Convert sample IDs from factor to character
        sample_ids <- as.character(
            DRIMSeq::samples(d)$sample_id
        )

        # Extract numeric expression matrix
        expr_gene <- as.matrix(
            gene_counts[, sample_ids, drop = FALSE]
        )

        # Calculate total gene expression per sample
        gene_expression <- colSums(expr_gene)

        # Calculate mean and median across all samples
        mean_expression <- mean(gene_expression)
        median_expression <- median(gene_expression)

        # Generate DRIMSeq plot
        plot <- DRIMSeq::plotProportions(
            x = d,
            gene_id = gene,
            group_variable = "condition",
            plot_type = "barplot",
            order_features = TRUE,
            order_samples = TRUE,
            plot_fit = TRUE,
            plot_main = TRUE
        )

        # Add expression statistics to title
        plot <- plot +
            ggtitle(
                paste0(
                    gene,
                    "\nMean expression = ",
                    round(mean_expression, 2),
                    ", Median expression = ",
                    round(median_expression, 2)
                )
            ) +
            theme_bw()

        # Save plot
        ggsave(
            filename = file.path(
                out_dir,
                paste0(
                    comparison,
                    "_gene_",
                    gene,
                    "_proportions.png"
                )
            ),
            plot = plot,
            width = 10,
            height = 6,
            dpi = 300
        )
    }

    # --------------------------------------------------
    # Combined overview plot:
    # Top transcripts ranked by adjusted p-value
    # --------------------------------------------------

    top_features$feature_label <- paste(
        top_features$gene_id,
        top_features$feature_id,
        sep = " | "
    )

    top_features$neg_log10_adj_p <- -log10(
        pmax(
            top_features$adj_pvalue,
            .Machine$double.xmin
        )
    )

    top_features$feature_label <- factor(
        top_features$feature_label,
        levels = rev(top_features$feature_label)
    )

    p <- ggplot(
        top_features,
        aes(
            x = neg_log10_adj_p,
            y = feature_label
        )
    ) +
        geom_col(fill = "steelblue") +
        labs(
            title = paste(
                "Top",
                top_n,
                "DTU transcript features:",
                comparison
            ),
            x = "-log10 adjusted p-value",
            y = "Gene ID | Transcript ID"
        ) +
        theme_bw() +
        theme(
            axis.text.y = element_text(size = 8),
            plot.title = element_text(
                hjust = 0.5
            )
        )

    ggsave(
        filename = file.path(
            out_dir,
            paste0(
                comparison,
                "_top",
                top_n,
                "_DTU_overview.png"
            )
        ),
        plot = p,
        width = 10,
        height = 7,
        dpi = 300
    )

    cat(
        "Plots saved for:",
        comparison,
        "\n"
    )

    return(top_features)
}

# ==================================================
# 4. Run plotting for all comparisons
# ==================================================

top10_results <- list()

for (comparison in names(results)) {

    top10_results[[comparison]] <- plot_top10_dtu(
        result_list = results[[comparison]],
        comparison = comparison,
        d = d,
        top_n = 10,
        padj_cutoff = 0.05
    )
}

# ==================================================
# 5. Combine top 10 results
# ==================================================

top10_all <- do.call(
    rbind,
    top10_results
)

if (!is.null(top10_all)) {

    write.csv(
        top10_all,
        file.path(
            out_dir,
            "top10_DTU_all_comparisons.csv"
        ),
        row.names = FALSE
    )
}

cat("\nAll DRIMSeq plots and results saved.\n")

