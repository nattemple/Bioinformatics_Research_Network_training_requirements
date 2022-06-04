## Load libraries
library("recount")
library(rlang)
library(stringr)
suppressMessages( library("DESeq2") )
suppressMessages( library(org.Hs.eg.db) )
suppressMessages( library( AnnotationDbi ) )
suppressMessages( library( dplyr ) )


# from https://jhubiostatistics.shinyapps.io/recount/
# 'IL-1ß and SERPINA3 are markers of an aggressive Barretts Oesophagus phenotype identified using mRNA sequencing'
accession='SRP043694'

## Download the RangedSummarizedExperiment object
url <- download_study(accession)

## View the url for the file by printing the object url
url

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

# temp
#for(i in 1:length(rse$characteristics)){
#  print(rse$characteristics[[i]][[1]])
#  #subject_ids[length(subject_ids)+1] <- as.numeric(substring(rse$characteristics[[i]][[1]],14,))
#}

# create factor cond: Simple intestinal metaplasia vs. Low Grade Dsyplasia vs. High Grade Dsyplasia
colData(rse) <- as.data.frame(colData(rse)) %>% 
  mutate(cond = substring(lapply(characteristics, `[[`, 1), 22)) %>%
  mutate(cond=str_replace_all(cond, " ", "_")) %>%
  DataFrame

# convert to factor
rse$cond <- factor(rse$cond)

# Differential Expression analysis
dds <- DESeqDataSet(rse, design = ~ cond)
dds <- DESeq(dds)
res <- results(dds, contrast = c("cond", "High_Grade_Dsyplasia", "Simple_intestinal_metaplasia"))

# sort by padj
res <- res[with(res, order(padj)), ]

# set up differentially expressed gene selection criteria 
pThr        <- 0.01   # for adj. p-value threshold (<= 0.01)
logFCThr    <- 1      # for log2 fold-change (at least 2 fold)
baseMeanThr <- 30     # for average expression level (at least 30 counts).
cpmThr      <- 1      # for copy-per-million (at least 1 cpm). 

# annotate sigRes
anno <- AnnotationDbi::select(org.Hs.eg.db, rownames( res ),
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
                              keytype="ENSEMBL")

res = cbind( ENSEMBL = rownames( res), res )
outTable <- left_join( as.data.frame( res ), anno )
head( outTable ) 

idx = which( outTable$padj <= pThr 
             & abs( outTable$log2FoldChange ) >= logFCThr 
             & outTable$baseMean >= baseMeanThr 
)
sigRes = outTable[idx, ]
sigRes

# PCA starts here

# perform regularized-logarithm transformation (rlog)
rld <- rlog(dds)

pca_plot <- plotPCA(rld, intgroup = c("cond"))
