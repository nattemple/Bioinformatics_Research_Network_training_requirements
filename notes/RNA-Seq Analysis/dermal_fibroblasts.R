## Load libraries
library(rlang)
library(stringr)
library(data.table)
library(DT)
library(dplyr)
suppressMessages(library("DESeq2"))
library(ggplot2)
suppressMessages(library(EnhancedVolcano))
library(pheatmap)
library(PoiClaClu)
library(RColorBrewer)
library(clusterProfiler)
suppressMessages(library("recount"))
library(EnsDb.Hsapiens.v86)
suppressMessages(library(AnnotationDbi))
library(msigdbr)

accession <- "SRP012461"

## Download the RangedSummarizedExperiment object
url <- download_study(accession)

## Load the data
load(file.path(accession, "rse_gene.Rdata"))

# truncate the ensemble id version number from the rownames
rownames <- rownames(rse_gene)
for (i in 1:length(rownames)) {
  rownames[i] <- gsub("\\..*", "", rownames[i])
}
rownames(rse_gene) <- rownames

# truncate also the gene_id column
rowdata <- rowData(rse_gene)
rowdata$gene_id <- gsub("\\..*", "", rowdata$gene_id)
rowData(rse_gene) <- rowdata

# assign condition from examining ColData(rse_gene)
rse_gene$cond <- c("untreated", "wt-myoD", "vp64-MyoD", "untreated", "wt-myoD", "vp64-MyoD")

# use raw counts
rse <- rse_gene

# convert to factor
rse$cond <- factor(rse$cond)

# Differential Expression analysis
dds <- DESeqDataSet(rse, design = ~cond)
dds <- DESeq(dds)

res <- results(dds, contrast = c("cond"
                                 , "vp64-MyoD"
                                 , "wt-myoD")) 

# sort by padj
res <- res[with(res, order(padj)), ]

# set up differentially expressed gene selection criteria
pThr <- 0.01 # for adj. p-value threshold (<= 0.01)
logFCThr <- 1 # for log2 fold-change (at least 2 fold)
baseMeanThr <- 30 # for average expression level (at least 30 counts).

# annotate sigRes (org.Hs.eg.db)
anno <- ensembldb::select(
  EnsDb.Hsapiens.v86,
  filter = GeneIdFilter("ENS", "startsWith"),
  keys = rownames(res),
  keytype = "GENEID",
  columns = c("SYMBOL", "GENEID", "GENENAME")
)

res <- cbind(ENSEMBL = rownames(res), res)
outTable <- left_join(as.data.frame(res), anno, by = c("ENSEMBL" = "GENEID"))


# significant genes
idx <- which(outTable$padj <= pThr &
               abs(outTable$log2FoldChange) >= logFCThr &
               outTable$baseMean >= baseMeanThr)
sigRes <- outTable[idx, ]

# significant genes (upregulated)
idx <- which(outTable$padj <= pThr &
               abs(outTable$log2FoldChange) >= logFCThr &
               outTable$log2FoldChange >= 0 &
               outTable$baseMean >= baseMeanThr)
upregulated <- outTable[idx, ]

# significant genes (downregulated)
idx <- which(outTable$padj <= pThr &
               abs(outTable$log2FoldChange) >= logFCThr &
               outTable$log2FoldChange <= 0 &
               outTable$baseMean >= baseMeanThr)
downregulated <- outTable[idx, ]

# GSEA

# Get the over-expressed genes as a vector
over_expressed_genes <- upregulated %>%
  pull(SYMBOL)

# Get the gene sets and wrangle
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5")
gene_sets <- gene_sets %>%
  dplyr::select(gs_name, gene_symbol)

# Run over-representation analysis
egmt <- enricher(
  gene = over_expressed_genes,
  TERM2GENE = gene_sets
)
edf <- as.data.frame(egmt)

# perform regularized-logarithm transformation (rlog)
# rld <- rlog(dds)
rld <- vst(dds)

save(sigRes, file = "sigRes.RData")
save(upregulated, file = "upregulated.RData")
save(downregulated, file = "downregulated.RData")
save(dds, file = "dds.RData")
save(rld, file = "rld.RData")
save(egmt, file = "egmt.RData")

# analyze PAX6 in wt-MyoD over untreated compared with vp64-MyoD over untreated

# wt-MyoD over untreated

res <- results(dds, contrast = c("cond"
                                 , "wt-myoD"
                                 , "untreated")) 

# sort by padj
res <- res[with(res, order(padj)), ]

# set up differentially expressed gene selection criteria 
pThr        <- 0.01   # for adj. p-value threshold (<= 0.01)
logFCThr    <- 1      # for log2 fold-change (at least 2 fold)
baseMeanThr <- 30     # for average expression level (at least 30 counts).

# annotate sigRes (org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
anno <- ensembldb::select(
  EnsDb.Hsapiens.v86, 
  filter = GeneIdFilter("ENS", "startsWith"),
  keys = rownames( res ), 
  keytype = "GENEID", 
  columns = c("SYMBOL","GENEID","GENENAME"))

res = cbind( ENSEMBL = rownames( res), res )
outTable <- left_join( as.data.frame( res ), anno, by = c("ENSEMBL" = "GENEID") )

wt_MyoD_over_untreated <-
  dplyr::filter(outTable, grepl("PAX",GENENAME) | grepl("MYOG",GENENAME))
# sort 
wt_MyoD_over_untreated <- wt_MyoD_over_untreated[with(wt_MyoD_over_untreated, order(GENENAME)), ]

save(egmt, file = "egmt.RData")
save(wt_MyoD_over_untreated, file = "wt_MyoD_over_untreated.RData")

# check hypothesis

# assign condition from examining ColData(rse_gene)
rse_gene$replicate <- c('1','1','1','2','2','2')

rse_gene$replicate <- factor(rse_gene$replicate)

dds <- DESeqDataSet(rse_gene, design = ~replicate)

dds <- DESeq(dds)

# check cell replicates as to which one has the higher PAX3

normData <- counts(dds, normalized=TRUE)

# annotate sigRes (org.Hs.eg.db)
anno <- ensembldb::select(
  EnsDb.Hsapiens.v86,
  filter = GeneIdFilter("ENS", "startsWith"),
  keys = rownames(normData),
  keytype = "GENEID",
  columns = c("SYMBOL", "GENEID", "GENENAME")
)

normData <- cbind(ENSEMBL = rownames(normData), normData)
normData <- left_join(as.data.frame(normData), anno, by = c("ENSEMBL" = "GENEID"))

normData <-
  dplyr::filter(normData, GENENAME=="PAX3")

# sort 
normData <- normData[with(normData, order(GENENAME)), ]

normData <- subset(normData,select=c("SRR1614280","SRR1614283", "GENENAME"))


# second replicate has higher initial PAX3
#SRR1614280       SRR1614283 GENENAME
#1 15617.6538464749 14864.7579475358     PAX3



res <- results(dds, contrast = c("replicate"
                                 , "2"
                                 , "1")) 

# annotate sigRes (org.Hs.eg.db)
anno <- ensembldb::select(
  EnsDb.Hsapiens.v86,
  filter = GeneIdFilter("ENS", "startsWith"),
  keys = rownames(res),
  keytype = "GENEID",
  columns = c("SYMBOL", "GENEID", "GENENAME")
)

res <- cbind(ENSEMBL = rownames(res), res)
outTable <- left_join(as.data.frame(res), anno, by = c("ENSEMBL" = "GENEID"))


outTable <-
  dplyr::filter(outTable, grepl("PAX",GENENAME) | grepl("MYOG",GENENAME))
# sort 
outTable <- outTable[with(outTable, order(GENENAME)), ]
