## Load libraries
library("recount")
library(rlang)
suppressMessages( library("DESeq2") )
suppressMessages( library(org.Hs.eg.db) )
suppressMessages( library( AnnotationDbi ) )
suppressMessages( library( dplyr ) )


# from https://jhubiostatistics.shinyapps.io/recount/
accession='SRP055514'

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

#temp
# try to annotate gene names directly from input data rse
#anno <- AnnotationDbi::select(org.Hs.eg.db, rownames( rse ),
#                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
#                              keytype="ENSEMBL")
#
#rse = cbind( ENSEMBL = rownames( rse), rse )
#outTable <- left_join( as.data.frame( rse ), anno )
#head( outTable ) 
# end - temp



# create factor time_point: T0 or T3 (before bariatric surgery vs. 3 months after bariatric surgery)
colData(rse) <- as.data.frame(colData(rse)) %>% 
  mutate(time_point = substring(lapply(characteristics, `[[`, 2), 13, 14)) %>%
  DataFrame

# create factor type_of_surgery: adjustable_gastric_banding vs. Roux-en-Y_gastric_bypass
colData(rse) <- as.data.frame(colData(rse)) %>% 
  mutate(type_of_surgery = substring(lapply(characteristics, `[[`, 3), 29)) %>%
  DataFrame

# convert to factors
rse$time_point <- factor(rse$time_point)
rse$type_of_surgery <- factor(rse$type_of_surgery)

# Differential Expression analysis
dds <- DESeqDataSet(rse, design = ~ time_point)
dds <- DESeq(dds)
res <- results(dds, contrast = c("time_point", "T3", "T0"))

# sort by padj
res <- res[with(res, order(padj)), ]

# set up differentially expressed gene selection criteria 
pThr        <- 0.01   # for adj. p-value threshold (<= 0.01)
logFCThr    <- .1      # for log2 fold-change (at least 2 fold)
baseMeanThr <- 20     # for average expression level (at least 20 counts).
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

