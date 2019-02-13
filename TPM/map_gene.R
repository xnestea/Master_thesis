setwd("Master/bin/")

library(biomaRt)
load("../../Data/TCGA data/uni_transcript_counts_exp_clinical_add1.Rdata")

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
### for all ensemble genes in dataframe add information about transcript biotype and gene symbol
map_gene = getBM(attributes=c('ensembl_transcript_id','ensembl_gene_id',"transcript_biotype",'external_gene_name'), 
                 filters = 'ensembl_transcript_id', values=rownames(exp), mart=ensembl)
### select only protein coding.
#save_mp <- map_gene
#map_gene <- save_mp
map_gene = map_gene[which(map_gene$transcript_biotype=="protein_coding"),]
### change name of column
colnames(map_gene)[colnames(map_gene)=="external_gene_name"] <- "Gene"
### set correct rownames
rownames(map_gene) <- map_gene[,1]
### remove redundand columns
map_gene[,c(1:3)] <- NULL
save(map_gene, file = "../../bin/TPM/map_gene.Rdata")