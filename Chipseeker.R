if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("Biostrings")


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

BiocManager::install("ReactomePA")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", force = TRUE)

install.packages("dplyr")

BiocManager::install("org.Dm.eg.db")
BiocManager::install("org.Hs.eg.db")


#load packages
library(ChIPseeker)
library(clusterProfiler)
library(Biostrings)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ReactomePA)
library(org.Dm.eg.db)
library("ggplot2")
library("ggrepel")
library("ggplot2")
library("ggrepel")
library("gridExtra")
library("tidyr")
library("reshape")
library("data.table")
library(dplyr)
library(openxlsx)

rm(list = ls())

H3 <- read.table(file= "GSM409308_UCSD.H1.H3K4me3.LL227.bed", header = FALSE)
names(H3)[names(H3) == "V1"] <- "chr"
names(H3)[names(H3) == "V2"] <- "start"
names(H3)[names(H3) == "V3"] <- "end"
names(H3)[names(H3) == "V4"] <- "score"
H3 <- select(H3, chr, start, end)
#H3$chr <- sub("^", "chr", H3$chr)
H3 = GenomicRanges::makeGRangesFromDataFrame(H3)

head(H3)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#heatmap of binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)

tagMatrix <- getTagMatrix(H3, windows=promoter)

tagHeatmap(tagMatrix, xlim=c(-1000, 1000), color="black")
ggsave("Test.pdf")

peakAnno <- annotatePeak(H3, tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnno)
ggsave("Test2.pdf")


peakAnnoList <- annotatePeak(H3, TxDb=txdb, 
                             tssRegion=c(-1000, 1000), verbose=FALSE)
peakAnnoList

H3_annot <- as.data.frame(peakAnnoList@anno)
write.xlsx(H3_annot, file="Test3.xlsx")
