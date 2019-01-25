library(DESeq2)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(biomaRt)
library(piano)
library(stringr)

setwd("C:/Users/xnestea/Documents/Master")
### download TCGA data
query <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Transcriptome Profiling",
                  legacy = FALSE,
                  workflow.type = "HTSeq - Counts"
)

GDCdownload(query)
data <- GDCprepare(query)

### save barcodes related to given groups for later use in DE 
ccA_sample <- data$barcode[which(data$subtype_mRNA_cluster == 1)]
ccB_sample <- data$barcode[which(data$subtype_mRNA_cluster %in% c(2,3))]

### formating
a <- data.frame("barcode" = ccA_sample, "type" = c("ccA"))
b <- data.frame("barcode" = ccB_sample, "type" = c("ccB"))
DS <- rbind(a, b)

### get TCGA data but in reasonable format I can actualy work on
df <- GDCprepare(query, 
                 save=TRUE,
                 save.filename = "expres.rda",
                 summarizedExperiment = FALSE)

### make rownames and drop first 5 lines -> data about the quality
df2 <- df[,-1]
rownames(df2) <- df[,1]
df2 <- df2[-(1:5), , drop = FALSE]

### select samples we are interested in
df2 <- subset(df2, select = c(ccA_sample, ccB_sample))

### as matrix to get numeric values for DE
counts <- as.matrix(df2)

### perform DESeq and write the output into file
deseq <- DESeqDataSetFromMatrix(countData=counts, colData=DS, design= ~type)
deseq <- DESeq(deseq)
res <- results(deseq)
df <- as.data.frame(res)
df <- df[complete.cases(df), ]
DE_down_reg <- df[df$log2FoldChange < 0,]
DE_down_reg <- DE_down_reg[DE_down_reg$padj < 0.05,]

DE_up_reg <- df[df$log2FoldChange > 0,]
DE_up_reg <- DE_up_reg[DE_up_reg$padj < 0.05,]
rownames(df) <- str_replace(rownames(df),
                                   pattern = ".[0-9]+$",
                                   replacement = "")
rownames(DE_up_reg) <- str_replace(rownames(DE_up_reg),
                            pattern = ".[0-9]+$",
                            replacement = "")
rownames(DE_down_reg) <- str_replace(rownames(DE_down_reg),
                            pattern = ".[0-9]+$",
                            replacement = "")

rownames(DE_up_reg)

write.table(df, file = "C:/Users/xnestea/Documents/Master/Data/deseq_TCGA_ccAccB.txt", sep = "\t")
write.table(DE_down_reg, file = "C:/Users/xnestea/Documents/Master/Data/deseq_TCGA_DE_down_reg.txt", sep = "\t")
write.table(DE_up_reg , file = "C:/Users/xnestea/Documents/Master/Data/deseq_TCGA_DE_up_reg.txt", sep = "\t")
### map ensembl ids to gene sy
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

prot_id = getBM(attributes=c('ensembl_gene_id',"transcript_biotype",'external_gene_name'), filters = "ensembl_gene_id", values=rownames(df), mart=ensembl)
prot_id = prot_id[which(prot_id$transcript_biotype=="protein_coding"),]
map_gene <- prot_id
map_gene <- getBM(attributes=c('ensembl_gene_id',
                               'external_gene_name'),
                  mart = ensembl)
colnames(map_gene)[colnames(map_gene)=="external_gene_name"] <- "Gene"
map_gene$Gene <- toupper(map_gene$Gene)
rownames(map_gene) <- map_gene[,1]
map_gene[,1] <- NULL
map_gene[,1] <- NULL


### genesets downloaded from go2msig
gset=loadGSC("Data/Homo_sapiens_GSEA_GO_sets_bp_symbols_April_2015.gmt")
DE=df[ ,c('log2FoldChange','pvalue')]
rownames(DE) <- str_replace(rownames(DE),
                      pattern = ".[0-9]+$",
                      replacement = "")
### make a new row in DE table for respective Gene symbol 
DE[,"Gene"] <- map_gene[row.names(DE), "Gene"]

### remove rows for which no gene symbols is available
DE <- DE[!is.na(DE$Gene),]

### remove rows with same gene symbols
DE <- DE[!duplicated(DE$Gene),]

###Gene symbols are now row names, replacing Ensembl entries
row.names(DE) <- DE$Gene

### we only need logfoldchange and pvalue as piano input
DE <- DE[ ,c('log2FoldChange','pvalue')] 

### make both logdfoldchange and p values as separate matrices
### with Genesymbols and rownames
pval <- as.matrix(DE[ ,2])
lfc <- as.matrix(DE[ ,1])
row.names(pval) <- row.names(DE)
row.names(lfc) <- row.names(DE)
### run GSEA and save results to file
gsaRes <- runGSA(pval,lfc,gsc=gset, nPerm = 10000, adjMethod = "fdr")
GSAsummaryTable(gsaRes, save = TRUE, file = "Data/TCGA_GSEA_ccA_ccB_all.txt")

is.recursive(gsaRes)

GSAheatmap(gsaRes, cutoff=25, adjusted=TRUE, ncharLabel=75, cellnote="none", columnnames="full", colorkey=TRUE, colorgrad=NULL, cex=0.4)
