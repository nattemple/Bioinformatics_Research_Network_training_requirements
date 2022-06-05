## Load libraries
library("recount")
library(rlang)
suppressMessages( library("DESeq2") )
suppressMessages( library(org.Hs.eg.db) )
suppressMessages( library( AnnotationDbi ) )
suppressMessages( library( dplyr ) )


# from https://jhubiostatistics.shinyapps.io/recount/
# 'Characterization of a network of tumor suppressor microRNA''s in T Cell acute lymphoblastic leukemia'
accession='SRP050223'

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

# create factor cond: 'tissue: normal thymus' vs. 'tissue: T cell acute lymphoblastic leukemia'
colData(rse) <- as.data.frame(colData(rse)) %>% 
  mutate(cond = substring(lapply(characteristics, `[[`, 1), 9, 14)) %>%
  #mutate(cond=str_replace(cond, " ", "_")) %>%
  DataFrame

# convert to factor
rse$cond <- factor(rse$cond)

# Differential Expression analysis
dds <- DESeqDataSet(rse, design = ~ cond)
dds <- DESeq(dds)
res <- results(dds, contrast = c("cond", "T cell", "normal"))

# sort by padj
res <- res[with(res, order(padj)), ]

# set up differentially expressed gene selection criteria 
pThr        <- 0.01   # for adj. p-value threshold (<= 0.01)
logFCThr    <- 4      # for log2 fold-change (at least 2 fold)
baseMeanThr <- 200     # for average expression level (at least 20 counts).
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


