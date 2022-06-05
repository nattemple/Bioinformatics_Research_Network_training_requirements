# from https://jhubiostatistics.shinyapps.io/recount/

## Load library
library("recount")


accession='SRP006676'
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

## Scale counts
#rse <- scale_counts(rse_gene)

# use raw counts
rse <- rse_gene

colData(rse)
rowData(rse)
head(assay(rse))


# create "cond" column with the 4 biological replicates:
# Healthy_non_smoker, Healthy_smoker, Smoker_with_Lung_Cancer, Smoker_without_Lung_Cancer
# pooling technical replicates Illumina and NuGEN (aka protocols)
library(dplyr)

rse$cancer <- case_when(
  startsWith(rse$title, "Healthy") ~ "No_cancer",
  grepl("without Lung Cancer", rse$title) ~ "No_cancer",
  TRUE ~ "Cancer"
)

rse$smoker <- case_when(
  grepl( "non-smoker", rse$title, fixed = TRUE) ~ "Non_smoker",
  TRUE ~ "Smoker"
)

rse$cond <- factor(paste0(rse$smoker, rse$cancer))
rse$smoker <- factor(rse$smoker)
rse$cancer <- factor(rse$cancer)

# create "protocol" column with the 2 protocols: Illumina, NuGEN:
rse$protocol <- case_when(
  grepl("Illumina", rse$title) ~ "Illumina",
  grepl("NuGEN", rse$title) ~ "NuGEN"
)

# check cross-tabulations
table(rse$title)
table(rse$title, rse$cond)
table(rse$cond)
table(rse$title, rse$cancer)
table(rse$title, rse$smoker)

suppressMessages( library("DESeq2") )
dds <- DESeqDataSet(rse, design = ~ cond)

dds <- DESeq(dds)

res <- results(dds)
res <- results(dds, contrast=c("cond", "SmokerNo_cancer", "SmokerCancer"))

#head(res)

# sort by padj
res <- res[with(res, order(padj)), ]

# We set up differentially expressed gene selection criterion here. 
pThr        <- 0.01   # for adj. p-value threshold (<= 0.01)
logFCThr    <- 1      # for log2 fold-change (at least 2 fold)
baseMeanThr <- 20     # for average expression level (at least 20 counts).
cpmThr      <- 1      # for copy-per-million (at least 1 cpm). 


idx = which( res$padj <= pThr & 
               abs( res$log2FoldChange ) >= logFCThr & 
               res$baseMean >= baseMeanThr )
sigRes = res[idx, ]
sigRes

suppressMessages( library(org.Hs.eg.db) )
suppressMessages( library( AnnotationDbi ) )
suppressMessages( library( dplyr ) )

#Check the column of human annotation library column names (ensembl ID?)
#columns(org.Hs.eg.db)

anno <- AnnotationDbi::select(org.Hs.eg.db, rownames( sigRes ),
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
                              keytype="ENSEMBL")

sigRes = cbind( ENSEMBL = rownames( sigRes), sigRes )
outTable <- left_join( as.data.frame( sigRes ), anno )
head( outTable ) 




# end



## Find a project of interest
project_info <- abstract_search("GSE14308")

project_info <- abstract_search("GSE32465")

## Download the gene-level RangedSummarizedExperiment data
download_study(project_info$project)


## Load the data
load(file.path(project_info$project, "rse_gene.Rdata"))

browseVignettes("DESeq2")


project_info <- abstract_search("GSE32465")
download_study(project_info$project)

# Airway
suppressMessages( library( airway ) )
data( "airway" )
se <- airway
head( assay( se ) )
rowData( se )
colData( se )
dds <- DESeqDataSet(se, design = ~ cell + dex)
dds <- DESeq(dds)
res <- results(dds)
pThr        <- 0.01   # for adj. p-value threshold (<= 0.01)
logFCThr    <- 1      # for log2 fold-change (at least 2 fold)
baseMeanThr <- 20     # for average expression level (at least 20 counts).
cpmThr      <- 1      # for copy-per-million (at least 1 cpm). 
idx = which( res$padj <= pThr & 
               abs( res$log2FoldChange ) >= logFCThr & 
               res$baseMean >= baseMeanThr )
sigRes = res[idx, ]
suppressMessages( library(org.Hs.eg.db) )
suppressMessages( library( AnnotationDbi ) )
suppressMessages( library( dplyr ) )

#Check the column of human annotation library column names (ensembl ID?)
columns(org.Hs.eg.db)
anno <- AnnotationDbi::select(org.Hs.eg.db, rownames( sigRes ), 
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
                              keytype="ENSEMBL")

sigRes = cbind( ENSEMBL = rownames( sigRes), sigRes )
outTable <- left_join( as.data.frame( sigRes ), anno )
head( outTable ) 


# test cbind
data_1 <- data.frame(x1 = c(7, 3, 2, 9, 0),          # Column 1 of data frame
                     x2 = c(4, 4, 1, 1, 8),          # Column 2 of data frame
                     x3 = c(5, 3, 9, 2, 4))          # Column 3 of data frame
y1 <- c(9, 8, 7, 6, 5) 
data_new1 <- cbind(data_1, y1)                       # cbind vector to data frame
data_new1 

# parathyroidSE
# https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf
BiocManager::install('parathyroidSE')
library( "parathyroidSE" )
data( "parathyroidGenesSE" )
dds <- DESeqDataSet(parathyroidGenesSE, design = ~ patient + treatment)

# alpineData
BiocManager::install('alpineData')
suppressMessages( library( alpineData ) )
data( "alpineData" )
colData(alpineData)
dds <- DESeqDataSet(alpineData, design = ~ patient + treatment)

# fission
BiocManager::install('fission')
suppressMessages( library( fission ) )
data( "fission" )
colData(fission)
dds <- DESeqDataSet(fission, design = ~ replicate + strain)

# macrophage
BiocManager::install('macrophage')
suppressMessages( library( macrophage ) )
data( "macrophage" )
colData(macrophage)
#dds <- DESeqDataSet(alpineData, design = ~ patient + treatment)

# oct4
BiocManager::install('oct4')
suppressMessages( library( oct4 ) )
dir <- system.file("extdata", package="oct4")
countdata <- assay( oct4 )
coldata <- read.csv(file.path(dir,"coldata.csv"))
dds <- DESeqDataSet(coldata, design = ~ patient + treatment)

# parathyroidSE
BiocManager::install('parathyroidSE')
suppressMessages( library( parathyroidSE ) )
countdata <- assay( parathyroidGenesSE )
coldata <- colData( parathyroidGenesSE )
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ patient + treatment)



library(GEOquery)
library(limma)
gse14308 <- getGEO("GSE14308", AnnotGPL = TRUE)[[1]]


gse <- getGEO("GSE33816", AnnotGPL = TRUE)[[1]]


# Load library
library("recount")

## Find a project of interest
project_info <- abstract_search("GSE32465")

download_study(project_info$project)

download_study(project_info$project)



library(recount)
project_info <- abstract_search('GSE66357')
download_study(project_info$project)
load(file.path(project_info$project, 'rse_gene.Rdata'))
rse_gene



## Install recount from Bioconductor
install.packages("BiocManager")
BiocManager::install('recount')

## Browse the vignetets for a quick description of how to use the package
library('recount')
browseVignettes('recount')

## Download the RangedSummarizedExperiment object at the gene level for 
## study SRP009615
url <- download_study('SRP009615')


