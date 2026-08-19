datadir = "/home/workscape/results/feature_counts"

library("sva")
library("ggplot2")
library("gridExtra")
library("edgeR")
library("UpSetR")
library("grid")


setwd(datadir)
uncorr_data = read.table("counts_all_feature.txt",header = TRUE, sep="\t")
head(uncorr_data)
dim(uncorr_data)

condition = c('B','C','C', 'A', 'A', 'B', 'C', 'B', 'A', 'A', 'A', 'A')
time = c('Y', 'Y', 'Y', 'X', 'Y', 'Y', 'Y', 'Y', 'Y', 'X', 'Y', 'X')
batch = c('pA', 'pC', 'pB', 'pA', 'pA', 'pB', 'pA', 'pC', 'pC', 'pB', 'pB', 'pC')
