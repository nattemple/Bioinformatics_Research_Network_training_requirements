## Load libraries
library(rlang)
library(stringr)
library(data.table)
library(DT)
library(dplyr)
suppressMessages( library("DESeq2") )
library(ggplot2)
suppressMessages( library( EnhancedVolcano ) )
library(pheatmap)
library(PoiClaClu)
library(RColorBrewer)

knitr::opts_chunk$set(echo = TRUE)

suppressMessages( library("recount") )
suppressMessages( library(org.Hs.eg.db) )
library(EnsDb.Hsapiens.v86)
suppressMessages( library( AnnotationDbi ) )

accession='SRP043694'

## Download the RangedSummarizedExperiment object
url <- download_study(accession)

## Load the data
load(file.path(accession, 'rse_gene.Rdata'))

# truncate the ensemble id version number from the rownames
rownames = rownames( rse_gene )
for(i in 1:length(rownames)){
  rownames[i] <- gsub("\\..*","",rownames[i])
}
rownames( rse_gene ) <- rownames

# truncate also the gene_id column
rowdata=rowData(rse_gene)
rowdata$gene_id = gsub("\\..*","",rowdata$gene_id)
rowData(rse_gene)=rowdata

# use raw counts
rse <- rse_gene


# Simple_intestinal_metaplasia
# High_Grade_Dysplasia
# Low_Grade_Dysplasia
colData(rse) <- as.data.frame(colData(rse)) %>% 
  mutate(cond = substring(lapply(characteristics, `[[`, 1), 22)) %>%
  mutate(cond=str_replace_all(cond, " ", "_")) %>%
  mutate(cond=str_replace_all(cond, "Dsyplasia", "Dysplasia")) %>% # correct typo
  DataFrame

# Barretts - focus on low vs. high dysplasia for now
rse <-
  subset(rse
         , select = colData(rse)$cond %in% 
           c("Simple_intestinal_metaplasia", "High_Grade_Dysplasia"))
colData(rse) <- as.data.frame(colData(rse)) %>% 
  # rename to Metaplasia vs. Dysplasia
  mutate(cond=str_replace_all(cond, "Simple_intestinal_metaplasia", "Metaplasia")) %>% 
  mutate(cond=str_replace_all(cond, "High_Grade_Dysplasia", "Dysplasia")) %>%
  DataFrame

# convert to factor
rse$cond <- factor(rse$cond)

# Differential Expression analysis
dds <- DESeqDataSet(rse, design = ~ cond)
dds <- DESeq(dds)

res <- results(dds, contrast = c("cond"
                                 , "Dysplasia"
                                 , "Metaplasia")) 

# sort by padj
res <- res[with(res, order(padj)), ]

# set up differentially expressed gene selection criteria 
pThr        <- 0.01   # for adj. p-value threshold (<= 0.01)
logFCThr    <- 1      # for log2 fold-change (at least 2 fold)
baseMeanThr <- 30     # for average expression level (at least 30 counts).

# annotate sigRes (org.Hs.eg.db)
#anno <- AnnotationDbi::select(org.Hs.eg.db, rownames( res ),
#                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
#                              keytype="ENSEMBL")
#res = cbind( ENSEMBL = rownames( res), res )
#outTable <- left_join( as.data.frame( res ), anno )

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


# significant genes
idx = which( outTable$padj <= pThr 
             & abs( outTable$log2FoldChange ) >= logFCThr 
             & outTable$baseMean >= baseMeanThr 
)
sigRes = outTable[idx, ]

# significant genes (upregulated)
idx = which( outTable$padj <= pThr 
             & abs( outTable$log2FoldChange ) >= logFCThr 
             & outTable$log2FoldChange >= 0 
             & outTable$baseMean >= baseMeanThr 
)
upregulated = outTable[idx, ]

# significant genes (downregulated)
idx = which( outTable$padj <= pThr 
             & abs( outTable$log2FoldChange ) >= logFCThr 
             & outTable$log2FoldChange <= 0 
             & outTable$baseMean >= baseMeanThr 
)
downregulated = outTable[idx, ]

write.csv(upregulated
          ,"\\Bioinformatics_Research_Network_training_requirements\\RNA-Seq Analysis\\upregulated.csv"
          , row.names = FALSE)

write.csv(downregulated
          ,"\\Bioinformatics_Research_Network_training_requirements\\RNA-Seq Analysis\\downregulated.csv"
          , row.names = FALSE)

# perform regularized-logarithm transformation (rlog)
#rld <- rlog(dds)
rld <- vst(dds)

# GSEA

library(msigdbr)
library(clusterProfiler)

# Get the over-expressed genes as a vector
 #optionally truncate version in SYMBOL 
 upregulated$SYMBOL <- gsub('\\..+', '', upregulated$SYMBOL, perl=T)
over_expressed_genes <- upregulated %>%
  pull(SYMBOL) 

# Get the gene sets and wrangle
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5")
gene_sets <- gene_sets %>%
  dplyr::select(gs_name, gene_symbol)

# Run over-representation analysis
egmt <- enricher(gene = over_expressed_genes,
                 TERM2GENE = gene_sets)
edf <- as.data.frame(egmt)

# Plot results with clusterProfiler
dotplot(egmt)
barplot(egmt)

