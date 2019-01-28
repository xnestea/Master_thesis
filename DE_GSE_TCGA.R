if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("piano", version = "3.8")

library(DESeq2)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(biomaRt)
library(piano)
library(stringr)


### set correct directory!!!!
setwd("C:/Users/xnestea/Documents/Master")


### RETRIEVE BARCODES FOR RESPECTIVE SUBSPECIES
### download TCGA data
query <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Transcriptome Profiling",
                  legacy = FALSE,
                  workflow.type = "HTSeq - Counts"
)

GDCdownload(query)
data <- GDCprepare(query)

### MAP TCGA SUBTYPES WITH CCA AND CCB
ccA_sample <- data$barcode[which(data$subtype_mRNA_cluster == 1)]
ccB_sample <- data$barcode[which(data$subtype_mRNA_cluster %in% c(2,3))]

### formating
a <- data.frame("barcode" = ccA_sample, "type" = c("ccA"))
b <- data.frame("barcode" = ccB_sample, "type" = c("ccB"))
### metadata file.
DS <- rbind(a, b)
rownames(DS) <- DS$barcode
DS[,1] <- NULL

### GET EXPRESSION DATA FROM TCGA AND FORMAT IT FOR DE
### get TCGA data but in reasonable format I can actualy work on
df <- GDCprepare(query, 
                 save=TRUE,
                 save.filename = "expres.rda",
                 summarizedExperiment = FALSE)

### make rownames and drop first 5 lines -> data about the quality
df2 <- df[,-1]
### set correct rownames and cut the endings of ensembl entries to be compatible with next step
rownames(df2) <- str_replace(df[,1],
                             pattern = ".[0-9]+$",
                             replacement = "")
df2 <- df2[-(1:5), , drop = FALSE]
### select samples we are interested in
df2 <- df2[, rownames(DS)]

### CONVERT ENSEMBLE IDS INTO GENE SUMBOLS AND SELECT ONLY PROTEIN CODING GENES
### Retreive human genes from ensembl
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
### for all ensemble genes in dataframe add information about transcript biotype and gene symbol
map_gene = getBM(attributes=c('ensembl_gene_id',"transcript_biotype",'external_gene_name'), filters = "ensembl_gene_id", values=rownames(df2), mart=ensembl)
### select only protein coding genes
map_gene = map_gene[which(map_gene$transcript_biotype=="protein_coding"),]
### change name of column
colnames(map_gene)[colnames(map_gene)=="external_gene_name"] <- "Gene"
#map_gene$Gene <- toupper(map_gene$Gene)
### set correct rownames
rownames(map_gene) <- map_gene[,1]
### remove redundand columns
map_gene[,c(1,2)] <- NULL
### make a new row in DE table for respective Gene symbol 
df3<- merge(df2, map_gene, by = "row.names", all.x=TRUE)
### remove rows for which no gene symbols is available
df3 <- df3[!is.na(df3$Gene),]
### remove rows with same gene symbols
df3 <- df3[!duplicated(df3$Gene),]
###Gene symbols are now row names, replacing Ensembl entries
row.names(df3) <- df3$Gene
### remove redundand columns
df3[,c("Gene","Row.names")] <- NULL
### as matrix to get numeric values for DE
df3 <- as.matrix(df3)

#write.table(df3, file = "C:/Users/xnestea/Documents/Master/Data/HT_seq_counts.txt", sep = "\t")


### DE ANALYSIS
### DE analysis
deseq <- DESeqDataSetFromMatrix(countData=df3, colData=DS, design= ~type)
deseq <- DESeq(deseq)
resultsNames(deseq)
res <- results(deseq)
### formating saving as data frame and removing rows containing NA
DE <- as.data.frame(res)
DE <- DE[complete.cases(DE), ]
### select significantly down regulated genes
DE_down_reg <- DE[DE$log2FoldChange < 0,]
DE_down_reg <- DE_down_reg[DE_down_reg$padj < 0.05,]
### select significantly up regulated genes
DE_up_reg <- DE[DE$log2FoldChange > 0,]
DE_up_reg <- DE_up_reg[DE_up_reg$padj < 0.05,]
### save files for further analysis
write.table(DE, file = "C:/Users/xnestea/Documents/Master/Data/deseq_TCGA_ccAccB_protcoding.txt", sep = "\t")
write.table(DE_down_reg, file = "C:/Users/xnestea/Documents/Master/Data/deseq_TCGA_DE_down_reg.txt", sep = "\t")
write.table(DE_up_reg , file = "C:/Users/xnestea/Documents/Master/Data/deseq_TCGA_DE_up_reg.txt", sep = "\t")

### GENE SET ENRICHMENT ANALYSIS (GSEA)
### genesets downloaded from go2msig
gset=loadGSC("Data/Homo_sapiens_GSEA_GO_sets_bp_symbols_April_2015.gmt")
### only those two values are necessary
DE<- DE[ ,c('log2FoldChange','pvalue')]

### make both logdfoldchange and p values as separate matrices
### with Genesymbols and rownames
pval <- as.matrix(DE[ ,2])
lfc <- as.matrix(DE[ ,1])
row.names(pval) <- row.names(DE)
row.names(lfc) <- row.names(DE)
### run GSEA and save results to file
gsaRes <- runGSA(pval,lfc,gsc=gset, nPerm = 10000, adjMethod = "fdr")
GSAsummaryTable(gsaRes, save = TRUE, file = "Data/TCGA_GSEA_ccA_ccB_all.txt")

### generate a heatmap depicting the results
GSAheatmap(gsaRes, cutoff=25, adjusted=TRUE, ncharLabel=75, cellnote="none", columnnames="full", colorkey=TRUE, colorgrad=NULL, cex=0.4)
